"""常染成分与监督判别模型。

设计原则（与 v05 一致）：
- 参考面板固定：目标个体不参与等位基因频率、位点选择、PCA、判别模型的拟合；
- 成分与判别均为"参考模型估计"，不是民族/国籍标签，也不是精确的遗传比例；
- 用交叉验证量化模型自身的分辨率（基线召回率、判别准确率），不乐观推断；
- 缺失不插补，仅使用完整调用位点。
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import nnls
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

SEED = 20261006


@dataclass
class Panel:
    """一个参考面板：包含哪些个体、用什么类别标签、MAF 阈值。"""

    name: str
    label_field: str
    categories: list[str]
    maf: float
    indices: np.ndarray
    labels: np.ndarray
    baseline: str = "em"  # 基线召回方式：em / nnls / none


@dataclass
class PanelResult:
    """一个面板对目标的完整估计结果。"""

    name: str
    categories: list[str]
    n_reference: int
    n_sites: int
    composition_em: dict[str, float]
    composition_nnls: dict[str, float]
    baseline_self_recall: dict[str, float]
    discriminator_accuracy: float
    discriminator_call: str
    discriminator_proba: dict[str, float]


def supervised_em(dosage: np.ndarray, freq: np.ndarray,
                  n_iter: int = 300, tol: float = 1e-9) -> np.ndarray:
    """监督 EM：固定各群体等位基因频率，估计单个个体的混合比例。

    freq: 类别 x 位点；dosage: 位点（0/1/2）。返回各类别比例（和为 1）。
    两个拷贝的来源独立，必须同时分配 alt 拷贝（g）与 ref 拷贝（2-g），
    否则中性高频变异的 alt 拷贝会被错误地平均稀释。
    """
    k = freq.shape[0]
    q = np.full(k, 1.0 / k)
    d = dosage[None, :]
    for _ in range(n_iter):
        # alt 拷贝（g 个）来自各类别的后验责任
        alt_num = q[:, None] * freq
        alt_den = alt_num.sum(axis=0)
        r_alt = np.divide(alt_num, alt_den, out=np.zeros_like(alt_num),
                          where=alt_den > 0)
        # ref 拷贝（2-g 个）来自各类别的后验责任
        ref_num = q[:, None] * (1.0 - freq)
        ref_den = ref_num.sum(axis=0)
        r_ref = np.divide(ref_num, ref_den, out=np.zeros_like(ref_num),
                          where=ref_den > 0)
        n_k = np.sum(r_alt * d, axis=1) + np.sum(r_ref * (2.0 - d), axis=1)
        total = n_k.sum()
        new_q = n_k / total if total > 0 else q
        if np.max(np.abs(new_q - q)) < tol:
            q = new_q
            break
        q = new_q
    return q


def nnls_mix(dosage: np.ndarray, freq: np.ndarray) -> np.ndarray:
    """非负最小二乘：dosage = 2 freq^T q，归一化。"""
    q, _ = nnls(2.0 * freq.T, dosage)
    s = q.sum()
    if s <= 0:
        return np.full(freq.shape[0], 1.0 / freq.shape[0])
    return q / s


def estimate_freq(genotypes: np.ndarray, labels: np.ndarray,
                  categories: list[str]) -> np.ndarray:
    """用参考个体估计各类别等位基因频率。"""
    freq = np.stack([genotypes[labels == c].mean(axis=0) / 2.0 for c in categories])
    return np.clip(freq, 1e-6, 1.0 - 1e-6)


def discriminator_cv(genotypes: np.ndarray, labels: np.ndarray,
                     target: np.ndarray, n_pc: int = 20):
    """PCA 降维 + LDA：返回 5 折准确率、目标判定与概率。"""
    model = make_pipeline(StandardScaler(), PCA(n_components=n_pc, random_state=SEED),
                          LinearDiscriminantAnalysis())
    cv = StratifiedKFold(5, shuffle=True, random_state=SEED)
    pred = cross_val_predict(model, genotypes, labels, cv=cv)
    accuracy = float((pred == labels).mean())
    model.fit(genotypes, labels)
    proba = model.predict_proba(target[None, :])[0]
    classes = list(model.classes_)
    call = classes[int(np.argmax(proba))]
    return accuracy, call, dict(zip(classes, map(float, proba)))


def baseline_self_recall(genotypes: np.ndarray, labels: np.ndarray,
                         categories: list[str], cap_per_class: int = 8) -> dict[str, float]:
    """成分模型自身召回率（5 折 CV，无泄漏）。

    每折用训练集估计频率，对测试集个体做 EM，统计每个真实类别平均判回自身的比例。
    为控制计算量，每类别每折最多抽 cap_per_class 个测试个体。
    """
    rng = np.random.default_rng(SEED)
    cv = StratifiedKFold(5, shuffle=True, random_state=SEED)
    sums = {c: 0.0 for c in categories}
    counts = {c: 0 for c in categories}
    for train, test in cv.split(genotypes, labels):
        freq = estimate_freq(genotypes[train], labels[train], categories)
        for c in categories:
            pool = test[labels[test] == c]
            if len(pool) == 0:
                continue
            chosen = pool if len(pool) <= cap_per_class else rng.choice(pool, cap_per_class, replace=False)
            for i in chosen:
                q = supervised_em(genotypes[i].astype(float), freq)
                sums[c] += float(q[categories.index(c)])
                counts[c] += 1
    return {c: (sums[c] / counts[c] if counts[c] else 0.0) for c in categories}


def baseline_self_recall_nnls(genotypes: np.ndarray, labels: np.ndarray,
                              categories: list[str]) -> dict[str, float]:
    """成分模型自身召回率（NNLS，5 折 CV，全量测试个体，速度快）。"""
    cv = StratifiedKFold(5, shuffle=True, random_state=SEED)
    sums = {c: 0.0 for c in categories}
    counts = {c: 0 for c in categories}
    for train, test in cv.split(genotypes, labels):
        freq = estimate_freq(genotypes[train], labels[train], categories)
        for i in test:
            q = nnls_mix(genotypes[i].astype(float), freq)
            c = labels[i]
            sums[c] += float(q[categories.index(c)])
            counts[c] += 1
    return {c: (sums[c] / counts[c] if counts[c] else 0.0) for c in categories}
