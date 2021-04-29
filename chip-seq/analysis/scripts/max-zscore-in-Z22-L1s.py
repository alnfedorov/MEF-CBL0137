import pathlib
from collections import defaultdict

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from pybedtools import BedTool
from tqdm import tqdm

import utils

# Gather L1Md A/T bound by Z22
blacklisted = BedTool(utils.paths.BLACKLIST)
z22 = BedTool(utils.paths.Z22_CURAX_14H).subtract(blacklisted, A=True)

repmasker = BedTool(utils.paths.REPEATMASKER)
repmasker = repmasker.subtract(blacklisted, A=True).intersect(z22, u=True, wa=True)

repeats = {"L1Md_A": [], "L1Md_T": []}
for repeat in tqdm(repmasker):
    if repeat.name in repeats:
        repeats[repeat.name].append(repeat)
repeats = {k: BedTool(v).sort() for k, v in repeats.items()}

# Get their sequences
sequences = {k: [] for k in repeats}
for repname, bed in repeats.items():
    fasta = bed.sort().sequence(
        fi=utils.paths.MM10_GENOME.as_posix(), s=True
    ).seqfn
    for seq in SeqIO.parse(fasta, format="fasta"):
        sequences[repname].append(str(seq.seq))


def maxz_segment(consensus, repeat, segmentation):
    matching = utils.miscellaneous.match_to_consensus(consensus, repeat)

    repeatz = utils.zhunt.run(repeat)
    repeatz = repeatz["Z-Score"].values.astype(np.float32)
    # determine a segment with the highest Z-score
    segments = {}

    # find orf1
    start, end = segmentation["ORF1"]
    orf1ind = matching[start: end]
    orf1ind = orf1ind[orf1ind >= 0]
    if len(orf1ind) > 0:
        segments['ORF1'] = repeatz[orf1ind].max()
        startorf1 = orf1ind[0]
        if startorf1 > 0:
            # 5`UTR is not truncated
            segments["5`UTR"] = repeatz[:startorf1].max()

    # find orf2
    start, end = segmentation["ORF2"]
    orf2ind = matching[start: end]
    orf2ind = orf2ind[orf2ind >= 0]
    if len(orf2ind) > 0:
        segments['ORF2'] = repeatz[orf2ind].max()
        endorf2 = orf2ind[-1]
        if endorf2 > 0:
            segments['3`UTR'] = repeatz[endorf2:].max()

    # ORF1 and ORF2 are absent -> either 5`UTR or 3`UTR
    if len(segments) == 0:
        where = np.flatnonzero(matching > 0)
        wmin, wmax, wmean = where.min(), where.max(), where.mean()
        z = repeatz[matching[matching > 0]].max()
        if wmean < segmentation["5`UTR"][1]:
            segments['5`UTR'] = z
        elif wmean > segmentation["3`UTR"][0]:
            segments['3`UTR'] = z
        else:
            assert False, "Unable to classify the repeat!"

    maxzseg = max(segments.items(), key=lambda x: x[1])
    return maxzseg


# Determine maxz segment for each repeat
maxz_segments = {"L1Md_A": [], "L1Md_T": []}
for repname, repeats in sequences.items():
    consensus = utils.l1consensus.sequence(repname)
    structure = utils.l1consensus.structure(repname)
    for sequence in repeats:
        maxz_segments[repname].append(maxz_segment(consensus, sequence, structure))

percentages = defaultdict(lambda: defaultdict(int))
for repname, segments in maxz_segments.items():
    for segment, _ in segments:
        percentages[repname][segment] += 1
    percentages[repname] = {k: v / len(segments) for k, v in percentages[repname].items()}

########################################################################################################################
# make the graph
colors = {
    "L1Md_T": "#F93B49",
    "L1Md_A": "#FF9554",
}
elemorder = ["L1Md_A", "L1Md_T"]
order = ["5`UTR", "ORF1", "ORF2", "3`UTR"]

barw = 0.5
gap = barw / 2
pad = barw / 2

anno_vertical_offset = 1.8
yticks = [15, 30, 45, 60, 75, 90]
lw = 0.75
fontsize = "large"

fig = plt.figure()
ax = fig.gca()
ax.set_xticks([])

ax.set_yticks(yticks)
ax.set_yticklabels([f"{y}%" for y in yticks], fontsize=fontsize)

offset = pad
for region in order:
    for ind, repname in enumerate(elemorder):
        value = round(percentages[repname][region] * 100, 1)
        rectangles = ax.bar(offset + ind * barw, value, color=colors[repname], width=barw,
                            edgecolor='black', align='edge')
        utils.miscellaneous.annotate_rectangles_with_values(rectangles, ax)
    ax.text(offset + len(colors) * barw / 2, -anno_vertical_offset, region,
            horizontalalignment='center', verticalalignment='top', fontsize=fontsize)
    offset = offset + len(colors) * barw + gap

ax.set_ylim((0, yticks[-1] * 1.1))
ax.set_xlim(0, offset)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)

legend = [mpatches.Patch(facecolor=colors[k], label=k.replace("_", " "),
                         edgecolor='black', linewidth=lw) for k in elemorder]
ax.legend(loc='upper right', handles=legend, frameon=False, fontsize=fontsize)
ax.set_title("Location of max z-score for L1Md A/L1Md T bound by Z22 (CBL0137 14h)",
             fontsize=fontsize).set_position([.5, 1.05])

saveto = pathlib.Path(__file__).name.replace(".py", ".svg")
fig.savefig(utils.paths.RESULTS.joinpath(saveto), bbox_inches="tight", pad_inches=0)
