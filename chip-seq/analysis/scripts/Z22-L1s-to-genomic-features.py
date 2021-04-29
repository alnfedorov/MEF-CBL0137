import os
import pathlib
import tempfile
from collections import defaultdict
from math import sin, radians, cos

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
from pybedtools import BedTool
from rpy2 import robjects

import utils


def annotate(bed: BedTool):
    fd, path = tempfile.mkstemp()
    os.close(fd)
    # Install and run R libraries
    # Not installed inside the docker to save time
    robjects.r(f"""
    # dirty hack to install packages on-demand
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    if (!require("ChIPseeker"))
        BiocManager::install("ChIPseeker")
    if (!require("TxDb.Mmusculus.UCSC.mm10.ensGene"))
        BiocManager::install("TxDb.Mmusculus.UCSC.mm10.ensGene")

    library(ChIPseeker)
    library(TxDb.Mmusculus.UCSC.mm10.ensGene)
    txdb = TxDb.Mmusculus.UCSC.mm10.ensGene
    peaks = annotatePeak("{bed.fn}", TxDb=txdb, tssRegion=c(-1000, 1000), verbose=TRUE)
    write.csv(peaks, file="{path}")
    """)
    df = pd.read_csv(path)
    os.remove(path)
    assert len(df) <= len(bed)

    # Merge some categories
    counts = df['annotation'].value_counts()
    finalcounts = defaultdict(int)
    for key, count in zip(counts.index, counts):
        if key.startswith("Intron") or key.startswith("Exon"):
            base = key.split(" ")[0]
            isfirst = key.count(" 1 of ")
            if isfirst:
                finalcounts[f"1st {base}"] += count
            else:
                finalcounts[base] += count
        else:
            finalcounts[key] += count
    # Normalization
    finalcounts = {k: v / sum(finalcounts.values()) for k, v in finalcounts.items()}
    return finalcounts


blacklisted = BedTool(utils.paths.BLACKLIST)
z22 = BedTool(utils.paths.Z22_CURAX_14H).subtract(blacklisted, A=True)

repeats = BedTool(utils.paths.REPEATMASKER)
repeats = repeats.subtract(blacklisted, A=True).intersect(z22, u=True, wa=True)

# gather and annotate target repeats
roi = {
    "L1Md_A": [],
    "L1Md_T": []
}
for r in repeats:
    if r.name in roi:
        roi[r.name].append(r)
roi = {k: annotate(BedTool(v)) for k, v in roi.items()}

# merge rare categories
for k, v in roi.items():
    for subk in ['1st Exon', "5' UTR", "3' UTR"]:
        if subk in v:
            v['Exon'] += v.pop(subk)
    for subk in ['Downstream (<1kb)', 'Downstream (2-3kb)', 'Downstream (1-2kb)']:
        if subk in v:
            v['Distal Intergenic'] += v.pop(subk)

# Set to 0 regions required categories if not present
order = ["Promoter", "Exon", "1st Intron", "Intron", "Distal Intergenic"]
for k, v in roi.items():
    for o in order:
        if o not in v:
            v[o] = 0
assert all(x in order for value in roi.values() for x in value)

########################################################################################################################
# make the plot
keys = [*order, *order[::-1]]
shares = [*(roi['L1Md_T'][v] for v in order), *(roi['L1Md_A'][v] for v in order[::-1])]

colorscheme = {
    "Distal Intergenic": "#DDDDDD",
    "Promoter": "#EEDD88",
    "Exon": "#EE8866",
    "1st Intron": "#99DDFF",
    "Intron": "#44BB99",
}

fig = plt.figure(figsize=(12, 8))
ax = fig.gca()

colors = [colorscheme[k] for k in keys]
wedges, _ = ax.pie(shares, startangle=-180,
                   wedgeprops=dict(width=0.5, edgecolor='#585A5B'), colors=colors)
ax.hlines(0, -1.3, 1.025, colors="black", linewidth=1)

ax.text(-1.2, 0.025, 'L1Md A', ha='center', va='bottom', fontsize='x-large')

angle = radians(180 - 35)
offset, offset2 = 1.1, 0.15
x, y = cos(angle) * offset, sin(angle) * offset
style = "Simple, tail_width=0.5, head_width=4, head_length=8"
kw = dict(arrowstyle=style, color='#585A5B')

ax.add_patch(
    mpatches.FancyArrowPatch((-offset, offset2), (x, y), connectionstyle=f"arc3,rad=-{offset / 10}", **kw)
)
ax.add_patch(
    mpatches.FancyArrowPatch((-offset, -offset2), (x, -y), connectionstyle=f"arc3,rad={offset / 10}", **kw)
)

ax.text(-1.2, -0.025, 'L1Md T', ha='center', va='top', fontsize='x-large')
ax.set_xlim(-1.5, 2)

legend = [mpatches.Patch(facecolor=colorscheme[k], label=k, edgecolor='black', linewidth=0.75) for k in order]
ax.legend(handles=legend, frameon=False, fontsize='x-large', loc='right')
ax.set_title("L1Md A/T bound by Z22 Ab\n(CBL0137-14h)", fontsize='xx-large', x=0, y=1.15, transform=ax.transData)
fig.tight_layout()

saveto = pathlib.Path(__file__).name.replace(".py", ".svg")
fig.savefig(utils.paths.RESULTS.joinpath(saveto), bbox_inches="tight", pad_inches=0)
