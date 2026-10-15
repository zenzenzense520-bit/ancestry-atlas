"""Single-file offline report. Every numeric result comes from computed artifacts."""
import base64
import csv
import html
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NAMES = {'AFR':'非洲','AMR':'美洲','EAS':'东亚','EUR':'欧洲','SAS':'南亚',
         'CHB':'北京汉族参考','CHS':'南方汉族参考','JPT':'东京日本人参考','CDX':'西双版纳傣族参考','KHV':'胡志明市京族参考'}
COLORS = {'AFR':'#786756','AMR':'#85939b','EAS':'#40584e','EUR':'#b0945e','SAS':'#b08276',
          'CHB':'#40584e','CHS':'#9c783c','JPT':'#85939b','CDX':'#b08276','KHV':'#786756'}
# 26 个 1000 Genomes 参考人群的中文说明（用于成分面板）。
NAMES.update({
    'ACB':'非洲加勒比（巴巴多斯）','ASW':'美国西南非裔','BEB':'孟加拉（达卡）','CEU':'犹他州西北欧裔',
    'CLM':'哥伦比亚（麦德林）','ESN':'尼日利亚埃桑','FIN':'芬兰人','GBR':'英国（英格兰/苏格兰）',
    'GIH':'古吉拉特印度裔（休斯敦）','GWD':'冈比亚','IBS':'伊比利亚（西班牙）','ITU':'印度泰卢固裔（英国）',
    'LWK':'肯尼亚卢希亚','MSL':'塞拉利昂门德','MXL':'墨西哥裔（洛杉矶）','PEL':'秘鲁（利马）',
    'PJL':'拉合尔旁遮普','PUR':'波多黎各','STU':'斯里兰卡泰米尔裔（英国）','TSI':'意大利托斯卡纳',
    'YRI':'尼日利亚约鲁巴'})

def esc(value): return html.escape(str(value))
def table(headers, rows, classname=''):
    return '<table class="'+classname+'"><thead><tr>'+''.join('<th scope="col">'+esc(x)+'</th>' for x in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+esc(x)+'</td>' for x in row)+'</tr>' for row in rows)+'</tbody></table>'

def inline_svg(path, prefix, classname, title):
    content=path.read_text(encoding='utf-8');content=content[content.index('<svg'):]
    for identifier in sorted(re.findall(r'\bid="([^"]+)"',content),key=len,reverse=True):
        content=content.replace('id="'+identifier+'"','id="'+prefix+identifier+'"')
        content=content.replace('#'+identifier+'"','#'+prefix+identifier+'"').replace('#'+identifier+')','#'+prefix+identifier+')')
    for group in NAMES:
        content=content.replace('id="'+prefix+'group-'+group+'"','id="'+prefix+'group-'+group+'" data-group="'+group+'"')
    content=content.replace('<svg ','<svg class="'+classname+'" role="img" aria-labelledby="'+prefix+'title" ',1)
    start=content.index('>')+1
    return content[:start]+'<title id="'+prefix+'title">'+esc(title)+'</title>'+content[start:]

def loco_text(p):
    return '、'.join(f'{g} {n}/{len(p["leave_one_chromosome_out"])} 次' for g,n in p['loco_nearest_counts'].items())

def projection_section(scope,p):
    east=scope=='east_asian';anchor='east' if east else 'projection';number='03' if east else '02'
    title='在东亚内部，再看一步。' if east else '把样本放回全球参考。'
    figures=''
    for pc in [2,3]:
        svgs=''.join(inline_svg(ROOT/'figures'/f'{scope}_pc1_pc{pc}{"_mobile" if mobile else ""}.svg',
            f'{scope}-pc{pc}{"m" if mobile else "d"}-','mobile-chart' if mobile else 'desktop-chart',
            f'{scope} 参考 PCA，PC1 对 PC{pc}。红色菱形为 Sample01。') for mobile in [False,True])
        figures+=f'<div class="chart-view" data-chart="{pc}" {"hidden" if pc==3 else ""}>{svgs}</div>'
    legend=''.join(f'<button type="button" data-group-toggle="{g}" aria-pressed="true" aria-label="显示{NAMES[g]}，{n}个体" style="--group-color:{COLORS[g]}"><i aria-hidden="true"></i>{g} <span>{n}</span></button>' for g,n in p['reference']['group_counts'].items())
    distances=table(['参考组','标准化距离'],[[f'{g} · {NAMES[g]}',f'{d:.3f}'] for g,d in sorted(p['center_distances'].items(),key=lambda x:x[1])],'distance-table')
    holdout=p['heldout_diagnostic'];weak=min(holdout['per_group'],key=lambda g:holdout['per_group'][g]['agreement'])
    diagnostic=table(['参考组','留出人数','一致率'],[[g,v['n'],f'{v["agreement"]*100:.2f}%'] for g,v in holdout['per_group'].items()])
    note=(f'<p class="east-note"><strong>细分结果的边界</strong>　{weak} 留出个体的最近中心一致率为 {holdout["per_group"][weak]["agreement"]*100:.1f}%。完整模型最近中心为 {p["nearest_reference_center"]}，80% 训练参考模型的样本最近中心为 {holdout["target_nearest_center"]}。相邻群中心排序对模型选择敏感，不能用来判定南北籍贯或具体身份。</p>' if east else '')
    variance=p['variance_explained'];n=p['reference']['reference_n']
    return f'''<section id="{anchor}" class="section" data-scope="{scope}" aria-labelledby="{anchor}-title"><div class="section-title"><div><div class="number">{number} / {"东亚内部投影" if east else "全球参考投影"}</div><h2 id="{anchor}-title">{title}</h2></div><span class="section-note">{n:,} 名参考 / {p['used_after_ld']:,} 个位点</span></div><p class="muted body-note">参考个体独立拟合，Sample01 仅投影。红色菱形为样本；坐标轴百分数表示解释的方差。{'CHB、CHS、JPT、CDX、KHV 为该发布中的五个东亚参考人群。' if east else 'AFR、AMR、EAS、EUR、SAS 为公开参考的五个粗分组。'}</p><div class="plot-controls"><div class="segmented" role="group" aria-label="{'东亚' if east else '全球'}主成分视图"><button type="button" data-view="2" aria-pressed="true">PC1 × PC2</button><button type="button" data-view="3" aria-pressed="false">PC1 × PC3</button></div><span class="view-caption small" aria-live="polite">PC1 {variance[0]*100:.2f}% · PC2 {variance[1]*100:.2f}% 解释方差</span></div><figure class="projection-figure">{figures}<figcaption>点击下方参考组切换显示。筛选只改变图形，计算使用完整对应参考。</figcaption></figure><div class="group-legend" aria-label="参考组显示开关">{legend}<span class="sample-legend"><i></i> Sample01</span></div>{note}<div class="grid2 projection-details"><div class="card"><div class="label">四维参考中心距离</div><h3>最近中心：{p['nearest_reference_center']}</h3><p class="small">PC1–PC4 按参考标准差缩放。距离是相对坐标比较，没有转成祖源比例或概率。两维最近中心为 {p['two_pc_nearest_center']}。</p>{distances}</div><div class="card soft"><div class="label">稳定性与适用范围</div><h3>{loco_text(p)}</h3><p class="small">逐次去掉一条常染色体后重新拟合。该检查描述模型敏感性，不能视为置信区间。</p><div class="diagnostic"><b>{holdout['centroid_agreement']*100:.2f}<small>%</small></b><span>参考内部留出集<br>最近中心与原分组的一致率</span></div><p class="small">80% 训练参考 / 20% 留出（{holdout['test_n']} 人），训练参考独立进行特征筛选。{weak} 一致率为 {holdout['per_group'][weak]['agreement']*100:.1f}%。这不是个人准确率。</p><details><summary>查看参考留出诊断</summary>{diagnostic}</details></div></div><div class="pipeline"><div><small>对齐候选位点</small><b>{p['usable_before_maf_ld']:,}</b></div><span aria-hidden="true">→</span><div><small>完整调用 / MAF 筛选</small><b>{p['after_maf']:,}</b></div><span aria-hidden="true">→</span><div><small>LD 筛选后</small><b>{p['used_after_ld']:,}</b></div><span aria-hidden="true">→</span><div><small>独立参考拟合</small><b>{n:,}</b></div></div></section>'''

COMP_TITLES={'global_super':'五个粗分组','global_populations':'26 个参考人群','east_asian':'东亚五个人群'}
COMP_DESC={'global_super':'AFR / AMR / EAS / EUR / SAS 五个粗分组，与全球 PCA 使用同一参考集合。',
 'global_populations':'在全部 26 个参考人群标签上重新估计，用于观察细分粒度下的稳定性。',
 'east_asian':'仅用 504 名东亚参考（CHB / CHS / JPT / CDX / KHV），与东亚 PCA 使用同一子集。'}

def comp_label(code):
    return f'{code} · {NAMES[code]}' if code in NAMES else str(code)

def comp_bars(panel,cutoff=.01,top=7):
    ordered=sorted(panel['composition_em'].items(),key=lambda kv:kv[1],reverse=True)
    head=[(k,v) for k,v in ordered if v>=cutoff][:top];taken={k for k,_ in head}
    tail=[(k,v) for k,v in ordered if k not in taken]
    rows=list(head)
    if tail:rows.append((f'其他 {len(tail)} 组',sum(v for _,v in tail)))
    return ''.join(f'<div class="comp-row"><span class="comp-name" title="{esc(comp_label(name))}">{esc(comp_label(name))}</span>'
                   f'<span class="comp-track"><i style="width:{value*100:.2f}%"></i></span>'
                   f'<span class="comp-val">{value*100:.1f}%</span></div>' for name,value in rows)

def comp_detail(panel):
    em=panel['composition_em'];nn=panel['composition_nnls']
    body=''.join(f'<tr><td>{esc(comp_label(k))}</td><td>{em[k]*100:.2f}%</td><td>{nn.get(k,0)*100:.2f}%</td></tr>'
                 for k,_ in sorted(em.items(),key=lambda kv:kv[1],reverse=True))
    return ('<details><summary>查看 EM 与 NNLS 逐项数值</summary><div class="table-scroll"><table class="comp-table">'
            '<thead><tr><th scope="col">类别</th><th scope="col">监督 EM</th><th scope="col">NNLS</th></tr></thead>'
            f'<tbody>{body}</tbody></table></div></details>')

def composition_section(c):
    cards=''
    for p in c['panels']:
        name=p['name'];call=p['discriminator_call'];recall=p['baseline_self_recall'].get(call,0.0)
        lo=min(p['baseline_self_recall'].values());hi=max(p['baseline_self_recall'].values())
        cards+=(f'<article class="card comp-panel"><div class="label">面板 · {esc(COMP_TITLES.get(name,name))}</div>'
                f'<h3>{esc(COMP_TITLES.get(name,name))}</h3>'
                f'<p class="small">{esc(COMP_DESC.get(name,""))} 参考 {p["n_reference"]:,} 人 / {p["n_sites"]:,} 位点。</p>'
                f'<p class="comp-legend small">条形为监督 EM 成分；悬停或展开可对比 NNLS。</p>'
                f'<div class="comp-bars">{comp_bars(p)}</div>'
                f'<div class="comp-foot">判别最近 <b>{esc(call)}</b> · 5 折 CV {p["discriminator_accuracy"]*100:.1f}% · '
                f'该类别自召回 {recall*100:.1f}%（面板 {lo*100:.0f}%–{hi*100:.0f}%）</div>{comp_detail(p)}</article>')
    lowest=min(min(p['baseline_self_recall'].values()) for p in c['panels'])
    note=(f'<p class="east-note"><strong>分辨率上限</strong>　成分比例来自固定参考频率的监督模型，不是民族或国籍标签。'
          f'参考个体自身的交叉验证召回率就是模型分辨率的上限：本报告面板中最低的类别自召回仅约 {lowest*100:.0f}%，'
          f'当某个已知参考个体都无法稳定判回自身时，不应把细分比例逐项当真。</p>')
    return (f'<section id="composition" class="section" aria-labelledby="composition-title">'
            f'<div class="section-title"><div><div class="number">04 / 常染成分与判别</div>'
            f'<h2 id="composition-title">在参考空间里，再估计一层成分。</h2></div>'
            f'<span class="section-note">固定参考频率 · 样本不参与拟合</span></div>'
            f'<p class="muted body-note">三个面板分别用固定参考等位基因频率的监督 EM（ADMIXTURE 风格）与非负最小二乘（NNLS）估计样本成分，'
            f'并用 PCA(20)+LDA 判别模型做 5 折交叉验证。等位基因频率、位点选择与 LD 均只由参考个体决定，Sample01 仅参与投影与估计。</p>'
            f'<div class="grid3 composition-panels">{cards}</div>{note}</section>')

def main():
    q=json.loads((ROOT/'results/sample_qc.json').read_text(encoding='utf-8'));b=json.loads((ROOT/'results/build_evidence.json').read_text(encoding='utf-8'))
    ref=json.loads((ROOT/'results/expanded_reference_qc.json').read_text(encoding='utf-8'))
    g=json.loads((ROOT/'results/global_pca.json').read_text(encoding='utf-8'));e=json.loads((ROOT/'results/east_asian_pca.json').read_text(encoding='utf-8'))
    c=json.loads((ROOT/'results/composition.json').read_text(encoding='utf-8'))
    fmt=lambda x:f'{x:,}'
    nearest=' 与 '.join(k for k,_ in sorted(e['center_distances'].items(),key=lambda x:x[1])[:2])
    states=q['coordinate_states'];raw=q['raw_call_kinds'];missing=raw.get('no_call',0)
    replacement={'INPUT_HASH':q['input_sha256'],'ALIGNED':fmt(ref['usable_before_maf_ld']),
       'EAST_NEIGHBORS':nearest,'GLOBAL_NEAREST':g['nearest_reference_center'],'GLOBAL_LOCO':'逐染色体复算为 '+loco_text(g),
       'RAW_ROWS':fmt(q['raw_rows']),'AUTOSOMAL':fmt(q['autosomal_clean_acgt_coordinates']),
       'GLOBAL_SITES':fmt(g['used_after_ld']),'EAST_SITES':fmt(e['used_after_ld']),
       'COORDINATES':fmt(q['unique_coordinates']),'DUP_GROUPS':fmt(q['coordinate_duplicate_groups']),
       'BASE_CONFLICTS':states.get('acgt_discordant',0),'MISSING_ROWS':fmt(missing),
       'MISSING_RATE':f'{missing/q["raw_rows"]*100:.2f}','ANCHOR_COUNT':len(b['anchors']),
       'COORD_HITS':fmt(ref['coordinate_hits']),'GLOBAL_SECTION':projection_section('global',g),'EAST_SECTION':projection_section('east_asian',e),
       'COMPOSITION_SECTION':composition_section(c)}
    for chrom,key in [('Y','Y_ROWS'),('MT','MT_ROWS')]:
        replacement[key]=next((r['raw_rows'] for r in q['chromosomes'] if r['chromosome']==chrom),0)
    replacement['Y_FINDING']='本次尚未进行 Y 单倍群判定。既有父系线索需要独立的 Y-SNP 或测序结果，当前无法核验下游分支。'
    replacement['QC_TABLE']=table(['指标','本文件结果'],[
        ['独立坐标',fmt(q['unique_coordinates'])],['重复坐标组',fmt(q['coordinate_duplicate_groups'])],
        ['重复坐标的额外行',fmt(q['redundant_coordinate_rows'])],['碱基共识坐标',fmt(q['clean_acgt_coordinates'])],
        ['常染碱基共识坐标',fmt(q['autosomal_clean_acgt_coordinates'])],['仅缺失坐标',fmt(states.get('no_call',0))],
        ['仅 I/D 编码坐标',fmt(states.get('indel_only',0))],['复杂编码排除组',len(q['duplicate_encoding_issues'])],
        ['A/C/G/T 冲突组',states.get('acgt_discordant',0)],['含后缀 rsID 行',fmt(q['suffixed_id_rows'])],['检测机构 / 芯片 / 链向','文件未注明']],'qc-summary')
    replacement['CHROMOSOME_TABLE']=table(['染色体','原始行','独立坐标','碱基共识','原始缺失','碱基共识杂合率'],[
        [r['chromosome'],fmt(r['raw_rows']),fmt(r['unique_coordinates']),fmt(r['clean_acgt_coordinates']),fmt(r['raw_no_call']),f'{r["clean_heterozygosity"]*100:.2f}%' if r['clean_heterozygosity'] is not None else '—'] for r in q['chromosomes']],'chromosome-table')
    replacement['ISSUE_TABLE']=table(['染色体','位置','源 ID','原始编码','排除原因'],[
        [r['chromosome'],fmt(r['position']),r['ids'].replace(';',' / '),r['source_genotypes'].replace(';',' / '),r['state']] for r in q['duplicate_encoding_issues']],'issue-table')
    replacement['ANCHOR_TABLE']=table(['位点','染色体','文件位置','GRCh37','GRCh38'],[[r['id'],r['chromosome'],fmt(r['observed_position']),fmt(r['grch37']),fmt(r['grch38'])] for r in b['anchors']],'anchor-table')
    labels={'direct_compatible':'直接兼容','unique_complement':'唯一反向互补','palindromic_excluded':'回文位点排除','no_unique_biallelic_snp':'没有唯一双等位 SNP','not_in_reference':'参考中未找到坐标','incompatible':'等位基因不兼容'}
    replacement['HARMONIZATION_TABLE']=table(['对齐状态','坐标数'],[[labels[k],fmt(v)] for k,v in ref['harmonization'].items()])
    replacement['SOFTWARE']=table(['软件','版本'],g['software'].items())
    replacement['REFERENCE_HASHES']=''.join('<p><b>'+esc(r['filename'])+'</b><br>'+esc(r['sha256'])+'</p>' for r in ref['download_manifest'])
    body=(ROOT/'src/report.body.html').read_text(encoding='utf-8')
    for token,value in replacement.items():body=body.replace('__'+token+'__',str(value))
    assert not re.search(r'__[A-Z_0-9]+__',body),'Unresolved placeholder'
    fonts=(ROOT/'src/fonts.css').read_text(encoding='utf-8')
    fonts=re.sub(r"url\('(fonts/[^']+)'\)",lambda m:"url('data:font/woff2;base64,"+base64.b64encode((ROOT/'src'/m[1]).read_bytes()).decode()+"')",fonts)
    css=fonts+'\n'+(ROOT/'src/base.css').read_text(encoding='utf-8')+'\n'+(ROOT/'src/report.css').read_text(encoding='utf-8')
    payload=json.dumps({'qc':q,'global':g,'east_asian':e,'composition':c,'assembly':b,'reference_qc':ref},ensure_ascii=False,separators=(',',':')).replace('<','\\u003c')
    notice='<!-- Font derived from Noto Serif SC; subset renamed Atlas Serif SC.\n'+(ROOT/'src/fonts/OFL.txt').read_text(encoding='utf-8').replace('--','—')+'\n-->'
    result='<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta name="color-scheme" content="light"><meta name="description" content="Sample01 完整 GRCh37 参考研究：全球与东亚独立 PCA、监督成分判别、调用质控及证据边界。"><title>Sample01 · 从全球走近东亚 · v0.5</title>'+notice+'<style>'+css+'</style></head><body>'+body+'<script type="application/json" id="sample-data">'+payload+'</script><script>'+(ROOT/'src/report.js').read_text(encoding='utf-8')+'</script></body></html>'
    (ROOT/'report/Sample01_reference_report_v05.html').write_text(result,encoding='utf-8')
    print(f'Built offline HTML: {len(result.encode()):,} bytes')

if __name__=='__main__':main()
