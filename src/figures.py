"""Publication/export figures from computed coordinates only; no invented points."""
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator

ROOT = Path(__file__).resolve().parents[1]
COLORS = {'AFR':'#786756', 'AMR':'#85939b', 'EAS':'#40584e', 'EUR':'#b0945e', 'SAS':'#b08276'}
COLORS.update({'CHB':'#40584e','CHS':'#9c783c','JPT':'#85939b','CDX':'#b08276','KHV':'#786756'})

def render(scope):
    result = json.loads((ROOT / f'results/{scope}_pca.json').read_text())
    rows = list(csv.DictReader((ROOT / f'results/{scope}_reference_pcs.tsv').open(), delimiter='\t'))
    coords = np.array([[float(r[f'PC{i}']) for i in range(1,5)] for r in rows])
    labels = np.array([r['group'] for r in rows])
    target = result['target_pcs']
    out = ROOT / 'figures'
    out.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':10, 'axes.labelcolor':'#6a675b',
                         'xtick.color':'#6a675b','ytick.color':'#6a675b','svg.fonttype':'path'})
    for second, mobile in [(1,False),(2,False),(1,True),(2,True)]:
        fig, ax = plt.subplots(figsize=(4.7,4.7) if mobile else (9,5.8), layout='constrained')
        fig.patch.set_facecolor('#f4efe4'); ax.set_facecolor('#f4efe4')
        for group in result['reference']['group_counts']:
            color = COLORS[group]
            chosen = labels == group
            artist = ax.scatter(coords[chosen,0], coords[chosen,second], s=11,
                                c=color, alpha=.48, edgecolors='none', label=f'{group} ({chosen.sum()})', zorder=2)
            artist.set_gid('group-' + group)
        sample = ax.scatter(target[0], target[second], s=145, c='#a54030', marker='D',
                            edgecolors='#f4efe4', linewidths=1.8, label='Sample01', zorder=5)
        sample.set_gid('sample-point')
        ax.annotate('Sample01', (target[0], target[second]), xytext=(13,-22), textcoords='offset points',
                    color='#a54030', weight='bold', fontsize=11, zorder=6,
                    bbox={'facecolor':'#f4efe4','edgecolor':'none','pad':3,'alpha':.9})
        ax.set_xlabel(f'PC1  |  {result["variance_explained"][0]*100:.2f}% variance', labelpad=12,fontsize=11 if mobile else 10)
        ax.set_ylabel(f'PC{second+1}  |  {result["variance_explained"][second]*100:.2f}% variance', labelpad=12,fontsize=11 if mobile else 10)
        ax.grid(color='#d6cebe', linewidth=.55, alpha=.7, zorder=1)
        for spine in ax.spines.values():spine.set_visible(False)
        ax.tick_params(length=0, pad=8,labelsize=11 if mobile else 10)
        if mobile:
            ax.xaxis.set_major_locator(MaxNLocator(4))
            ax.yaxis.set_major_locator(MaxNLocator(5))
        ax.margins(.12)
        name = f'{scope}_pc1_pc{second+1}' + ('_mobile' if mobile else '')
        fig.savefig(out / (name+'.svg'), metadata={'Date':None,'Creator':f'Reference-only PCA; {result["used_after_ld"]} sites'})
        # HTML provides interactive group controls; standalone PNG needs its own legend.
        handles, names = ax.get_legend_handles_labels()
        fig.legend(handles,names,loc='outside upper center',ncol=3 if mobile else 6,
                   frameon=False,fontsize=8 if mobile else 9,labelcolor='#6a675b')
        fig.savefig(out / (name+'.png'), dpi=180, metadata={'Software':'Sample01 reference-only PCA'})
        plt.close(fig)
    print(f'Rendered {scope} projections in SVG and PNG')

def main():
    for scope in ['global','east_asian']:
        render(scope)

if __name__ == '__main__':main()
