"""运行常染成分与监督判别模型：三个参考面板，输出 composition.json 与模型。"""
from __future__ import annotations

import argparse
import csv
import json
import platform
from pathlib import Path

import numpy as np
import sklearn

import composition as comp
from expanded_pca import select_sites

ROOT = Path(__file__).resolve().parents[1]
SEED = comp.SEED


def load_inputs():
    results = ROOT / "results"
    data = np.load(results / "reference_overlap.npz")
    samples = list(csv.DictReader((results / "reference_samples.tsv").open(), delimiter="\t"))
    sites = list(csv.DictReader((results / "reference_overlap_sites.tsv").open(), delimiter="\t"))
    g_all = data["genotypes"].T  # 个体 x 位点
    target = data["target_dosages"]
    return g_all, target, samples, sites, results


def build_panels(samples: list[dict]) -> list[comp.Panel]:
    n = len(samples)
    sup = np.array([r["group"] for r in samples])
    pop = np.array([r["population"] for r in samples])
    super_cats = ["AFR", "AMR", "EAS", "EUR", "SAS"]
    pop_cats = sorted(set(pop))
    eas_cats = ["CHB", "CHS", "JPT", "CDX", "KHV"]
    eas_idx = np.flatnonzero(np.isin(pop, eas_cats))
    return [
        comp.Panel("global_super", "group", super_cats, 0.01, np.arange(n), sup, baseline="em"),
        comp.Panel("global_populations", "population", pop_cats, 0.01, np.arange(n), pop,
                   baseline="nnls"),
        comp.Panel("east_asian", "population", eas_cats, 0.05, eas_idx, pop[eas_idx],
                   baseline="em"),
    ]


def run_panel(g_all, target, sites, panel: comp.Panel) -> tuple[comp.PanelResult, dict]:
    g_panel = g_all[panel.indices]
    selected, _ = select_sites(g_panel, sites, maf_min=panel.maf)
    if len(selected) < 100:
        raise ValueError(f"{panel.name}: 过滤后位点少于 100")
    g_sel = g_panel[:, selected]
    t_sel = target[selected]
    freq = comp.estimate_freq(g_sel, panel.labels, panel.categories)
    # 目标成分（EM、NNLS）
    q_em = comp.supervised_em(t_sel.astype(float), freq)
    q_nnls = comp.nnls_mix(t_sel.astype(float), freq)
    # 判别模型（PCA+LDA）
    accuracy, call, proba = comp.discriminator_cv(g_sel, panel.labels, t_sel)
    # 成分基线召回（交叉验证，无泄漏；按面板选择算法）
    if panel.baseline == "em":
        recall = comp.baseline_self_recall(g_sel, panel.labels, panel.categories)
    elif panel.baseline == "nnls":
        recall = comp.baseline_self_recall_nnls(g_sel, panel.labels, panel.categories)
    else:
        recall = {}
    print(f"{panel.name}: {len(selected)} 位点 / {len(g_panel)} 个体 完成", flush=True)
    result = comp.PanelResult(
        name=panel.name, categories=panel.categories, n_reference=len(g_panel),
        n_sites=int(len(selected)),
        composition_em=dict(zip(panel.categories, map(float, q_em))),
        composition_nnls=dict(zip(panel.categories, map(float, q_nnls))),
        baseline_self_recall=recall,
        discriminator_accuracy=accuracy, discriminator_call=call,
        discriminator_proba=proba)
    artifact = {"selected": selected, "freq": freq, "categories": panel.categories}
    return result, artifact


def build_output(panel_results: list[comp.PanelResult]) -> dict:
    return {
        "sample_alias": "Sample01", "version": "0.5",
        "reference": {"source": "PLINK public 1000 Genomes Phase 3 GRCh37",
                      "actual_ids_and_official_labels_verified": True},
        "method": {
            "fit": "Reference only; target never affects frequencies, site selection, PCA or discriminator",
            "composition": "Supervised EM (ADMIXTURE-style) and NNLS, reference allele frequencies fixed",
            "discriminator": "PCA(20) + LDA, 5-fold stratified cross-validation",
            "baseline": "Self-recall of known reference individuals via 5-fold CV (no leakage)",
            "ld_window_bp": 500000, "ld_r2_max": 0.2, "seed": SEED},
        "panels": [vars(pr) for pr in panel_results],
        "software": {"python": platform.python_version(), "numpy": np.__version__,
                     "sklearn": sklearn.__version__},
        "limitations": [
            "Composition is a reference-model estimate, not an ethnicity/nationality label.",
            "Ancestry proportions are sensitive to reference choice, LD and shared variation.",
            "Baseline self-recall is the model's resolution ceiling: do not read fine proportions literally "
            "when a known reference individual also does not return at 100%.",
            "Cross-validation accuracy describes the reference model, not a personal confidence interval.",
            "No ancient-DNA affinity, Y or mitochondrial haplogroup inference is included here."],
    }


def write_output(results: Path, panel_results: list[comp.PanelResult], artifacts: dict):
    """增量写盘：每个面板完成后立即写，避免中途中断导致全部丢失。"""
    (results / "composition.json").write_text(
        json.dumps(build_output(panel_results), indent=2, ensure_ascii=False), encoding="utf-8")
    np.savez_compressed(
        results / "composition_model.npz",
        **{f"{p}_selected": artifacts[p]["selected"] for p in artifacts},
        **{f"{p}_freq": artifacts[p]["freq"] for p in artifacts})


def run():
    g_all, target, samples, sites, results = load_inputs()
    panels = build_panels(samples)
    panel_results, artifacts = [], {}
    for panel in panels:
        pr, art = run_panel(g_all, target, sites, panel)
        panel_results.append(pr)
        artifacts[panel.name] = art
        write_output(results, panel_results, artifacts)  # 每个面板完成即写盘
        print(f"  -> 已增量写入 composition.json（{len(panel_results)}/{len(panels)}）", flush=True)
    for pr in panel_results:
        print(f"\n{pr.name}: 判别 {pr.discriminator_call} "
              f"(CV {pr.discriminator_accuracy:.3f})")
    print("\ncomposition.json 已写入", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.parse_args()
    run()
