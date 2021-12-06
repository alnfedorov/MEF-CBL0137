from pathlib import Path

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import utils

plt.rcParams['svg.fonttype'] = 'none'
sns.set_style("whitegrid")

PADJ = utils.DEG.RIP_PADJ
LFC = utils.DEG.RIP_LFC
TIMEPOINT = "48h"

df = pd.read_csv(utils.paths.DE_GENES.joinpath(f"{TIMEPOINT}_Z22_vs_{TIMEPOINT}_input.csv"))
df['-log10(padj)'] = -np.log10(df['padj'])
df = df[df['Biotype'] == 'protein_coding']

# Remove IgG enriched genes
IgG = utils.DEG.IgG()
df = df[~df['Ensembl Gene ID'].isin(IgG)]

# Mark ISG
ISG = utils.DEG.ISG('protein_coding')
df['ISG'] = df['Ensembl Gene ID'].isin(ISG)

# Mark edited genes
editing = utils.paths.RESULTS.joinpath("edited-genes.tsv")
edgenes = set(pd.read_csv(editing, sep='\t')['Gene'])
df['Edited'] = df['Gene name'].isin(edgenes)

# Make plot
LABELS_SIZE = 17
TITLE_SIZE = 20
TICKS_SIZE = 14

YLIM = 0, 32.5
XLIM = 0.95, 3.25

EDGEWIDTH = 1.5

fig = plt.figure(figsize=(12, 12))
ax = fig.gca()

subsets = [df[(~df.ISG) & (~df.Edited)],
           df[(df.ISG) & (~df.Edited)],
           df[(~df.ISG) & (df.Edited)],
           df[(df.ISG) & (df.Edited)]]
# Sanity check
genes = [set(x['Ensembl Gene ID']) for x in subsets]
assert all(len(genes[i].intersection(genes[j])) == 0 for i in range(len(genes)) for j in range(i + 1, len(genes)))

facecolors = ["#d7d7d7", "#F78969", "#d7d7d7", "#F78969"]
edgecolors = ["#b9b9b9", "#F74300", "#000000", "#000000"]
markers = ['o', 'o', 'o', 'o']
facealpha = [0.4, 0.4, 0.8, 0.8]
edgealpha = [0.65, 0.65, 1, 1]
for face, line, marker, facealpha, edgealpha, subdf in zip(facecolors, edgecolors, markers,
                                                           facealpha, edgealpha, subsets):
    ax.scatter(subdf['log2FoldChange'], subdf['-log10(padj)'], alpha=facealpha, linewidths=0, color=face,
               marker=marker, s=90)
    if edgealpha > 1e-3:
        ax.scatter(subdf['log2FoldChange'], subdf['-log10(padj)'], alpha=edgealpha, edgecolors=line,
                   linewidths=EDGEWIDTH, marker=marker, s=90).set_facecolor("none")

ax.tick_params(axis='both', which='major', labelsize=TICKS_SIZE)
ax.tick_params(axis='both', which='minor', labelsize=TICKS_SIZE)

# ax.set_ylabel("$-log_{10}(p$-$adj)$", fontsize=LABELS_SIZE)
ax.set_ylim(*YLIM)
ax.axhline(-np.log10(PADJ), color="black", linestyle="--")
# ax.set_xlabel("$log_2(fold~change)$", fontsize=LABELS_SIZE)
ax.set_xlim(*XLIM)
ax.axvline(LFC, color="black", linestyle="--")

for name in ['Eif2ak2', 'Ddx58', 'Ifih1', 'Xrn1', 'Slfn5', 'Agl', 'Knl1']:
    gene = df[df['Gene name'] == name]
    assert len(gene) == 1
    gene = gene.iloc[0]
    x, y = gene['log2FoldChange'], gene['-log10(padj)']
    ax.annotate(name, (x, y), xytext=(x + 0.75, y + 1.95), fontsize=LABELS_SIZE,
                arrowprops=dict(facecolor='#363636', edgecolor='#363636', arrowstyle='->', lw=2.5))

df = df[(df.padj < PADJ) & (df.log2FoldChange > LFC)]
labels = [
    f"non-ISG & not-edited ({((~df.ISG) & (~df.Edited)).sum()} Z22 DEGs)",
    f"ISG & not-edited ({((df.ISG) & (~df.Edited)).sum()} Z22 DEGs)",
    f"non-ISG & edited ({((~df.ISG) & (df.Edited)).sum()} Z22 DEGs)",
    f"ISG & edited ({((df.ISG) & (df.Edited)).sum()} Z22 DEGs)",
]
legend = [
    mlines.Line2D([], [], markerfacecolor=face, markeredgecolor=edge, marker=marker,
                  linestyle='None', markersize=TICKS_SIZE, label=lbl, markeredgewidth=EDGEWIDTH)
    for lbl, face, edge, marker in zip(labels, facecolors, edgecolors, markers)
]
legend = ax.legend(handles=legend, frameon=False, fontsize=TICKS_SIZE, loc='upper right', handletextpad=0.1)

x, y = ax.get_xlim()[1], ax.get_ylim()[1]
ax.text(x - 0.6, -np.log10(PADJ) * 1.1, f"p-adj < {PADJ}", fontsize=TICKS_SIZE, style='italic',
        weight='semibold')
ax.text(LFC * 1.02, y - 1, f"fold change > {2 ** LFC}", fontsize=TICKS_SIZE, style='italic', weight='semibold')

ax.set_title("$\\bf{Protein~coding~genes~differentially~enriched~in~Z22~RIP}$\n"
             f"(IFNβ-{TIMEPOINT}: Z22 vs Input, ADAR1-/- ZBP1-/-, {len(df)} Z22 DEGs)",
             fontsize=TITLE_SIZE)

print(TIMEPOINT)
print(f"\tISGs among Z22 DEGs: {df.ISG.sum() / len(df) * 100: .1f}%")
print(f"\tZ22 DEGs among ISGs: {df.ISG.sum() / len(ISG) * 100: .1f}%")

ax.xaxis.set_ticklabels(["" for _ in ax.xaxis.get_ticklabels()])
ax.yaxis.set_ticklabels(["" for _ in ax.yaxis.get_ticklabels()])

saveto = Path(__file__).name.replace(".py", f".png")
plt.savefig(utils.paths.RESULTS.joinpath(saveto), dpi=400, transparent=True, bbox_inches="tight", pad_inches=0)
