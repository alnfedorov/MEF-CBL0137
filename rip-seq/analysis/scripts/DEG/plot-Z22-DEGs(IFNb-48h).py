from pathlib import Path

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import utils

plt.rcParams['svg.fonttype'] = 'none'
sns.set_style("whitegrid")

LEGEND = True

PADJ = utils.DEG.RIP_PADJ
LFC = utils.DEG.RIP_LFC
TIMEPOINT = "48h"

df = pd.read_csv(utils.paths.DE_GENES.joinpath(f"{TIMEPOINT}_Z22_vs_{TIMEPOINT}_input.csv"))
df['-log10(padj)'] = -np.log10(df['padj'])
df = df[df['Biotype'] == 'protein_coding']

# Remove IgG enriched genes
IgG = utils.DEG.IgG()
df = df[~df['Ensembl Gene ID'].isin(IgG)]

# Select coding ISGs only
ISG = utils.DEG.ISG('protein_coding')
df = df[df['Ensembl Gene ID'].isin(ISG)]

# Mark edited genes
editing = utils.paths.RESULTS.joinpath("edited-genes.tsv")
editing = editing[editing['EditedIn'].isin({"utr3", "utr5", "exons"})]
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

subsets = [(~df.Edited), df.Edited]
subsets = [df[x] for x in subsets]

# Sanity check
genes = [set(x['Ensembl Gene ID']) for x in subsets]
assert all(len(genes[i].intersection(genes[j])) == 0 for i in range(len(genes)) for j in range(i + 1, len(genes)))

facecolors = ["#F78969", "#F78969"]
edgecolors = ["#F74300", "#000000"]
facealpha = [0.4, 0.8]
edgealpha = [0.4, 1]
for face, line, facealpha, edgealpha, subdf in zip(facecolors, edgecolors, facealpha, edgealpha, subsets):
    ax.scatter(subdf['log2FoldChange'], subdf['-log10(padj)'], alpha=facealpha, linewidths=0, color=face,
               marker='o', s=90)
    if edgealpha > 1e-3:
        ax.scatter(subdf['log2FoldChange'], subdf['-log10(padj)'], alpha=edgealpha, edgecolors=line,
                   linewidths=EDGEWIDTH, marker='o', s=90).set_facecolor("none")

ax.tick_params(axis='both', which='major', labelsize=TICKS_SIZE)
ax.tick_params(axis='both', which='minor', labelsize=TICKS_SIZE)

ax.set_ylim(*YLIM)
ax.axhline(-np.log10(PADJ), color="black", linestyle="--")
ax.set_xlim(*XLIM)
ax.axvline(LFC, color="black", linestyle="--")

if LEGEND:
    ax.set_ylabel("$-log_{10}(p$-$adj)$", fontsize=LABELS_SIZE)
    ax.set_xlabel("$log_2(fold~change)$", fontsize=LABELS_SIZE)

    for name in ['Eif2ak2', 'Ddx58', 'Ifih1', 'Xrn1', 'Slfn5', 'Knl1']:
        gene = df[df['Gene name'] == name]
        assert len(gene) == 1, name
        gene = gene.iloc[0]
        x, y = gene['log2FoldChange'], gene['-log10(padj)']
        ax.annotate(name, (x, y), xytext=(x + 0.75, y + 1.95), fontsize=LABELS_SIZE,
                    arrowprops=dict(facecolor='#363636', edgecolor='#363636', arrowstyle='->', lw=2.5))

    df = df[(df.padj < PADJ) & (df.log2FoldChange > LFC)]
    labels = [
        f"Not edited ({(~df.Edited).sum()} Z-forming ISG mRNAs)",
        f"Edited ({df.Edited.sum()} Z-forming ISG mRNAs)",
    ]
    legend = [
        mlines.Line2D([], [], markerfacecolor=face, markeredgecolor=edge, marker='o',
                      linestyle='None', markersize=TICKS_SIZE, label=lbl, markeredgewidth=EDGEWIDTH)
        for lbl, face, edge in zip(labels, facecolors, edgecolors)
    ]
    legend = ax.legend(handles=legend, frameon=False, fontsize=TICKS_SIZE, loc='upper right', handletextpad=0.1)

    x, y = ax.get_xlim()[1], ax.get_ylim()[1]
    ax.text(x - 0.6, -np.log10(PADJ) * 1.1, f"p-adj < {PADJ}", fontsize=TICKS_SIZE, style='italic',
            weight='semibold')
    ax.text(LFC * 1.02, y - 1, f"fold change > {2 ** LFC}", fontsize=TICKS_SIZE, style='italic', weight='semibold')

    ax.set_title("Coding ISGs differentially enriched in Z22 RIP-seq\n"
                 f"(IFNβ-{TIMEPOINT}: Z22 vs Input, ADAR1/ZBP1 dKO, {len(df)} Z22 DEGs)",
                 fontsize=TITLE_SIZE)

if not LEGEND:
    ax.xaxis.set_ticklabels(["" for _ in ax.xaxis.get_ticklabels()])
    ax.yaxis.set_ticklabels(["" for _ in ax.yaxis.get_ticklabels()])

saveto = Path(__file__).name.replace(".py", f"(legend={LEGEND}).png")
plt.savefig(utils.paths.RESULTS.joinpath(saveto), dpi=600, transparent=True, bbox_inches="tight", pad_inches=0)

saveto = Path(__file__).name.replace(".py", f"(legend={LEGEND}).tiff")
plt.savefig(utils.paths.RESULTS.joinpath(saveto), dpi=600, transparent=True, bbox_inches="tight", pad_inches=0)
