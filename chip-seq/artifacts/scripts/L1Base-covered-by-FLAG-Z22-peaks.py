import pathlib

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from pybedtools import BedTool

import utils

blacklisted = BedTool(utils.paths.BLACKLIST).sort()

orf1_orf2_intact = BedTool(utils.paths.L1BASE_ORF1_ORF2_INTACT).sort()
orf2_intact = BedTool(utils.paths.L1BASE_ORF2_INTACT).sort()
full_length = BedTool(utils.paths.L1BASE_FULL_LEN_NON_INTACT).sort()

# remove intersections
l1s = {
    "ORF1 & ORF2 intact": orf1_orf2_intact,
    "ORF2 intact": orf2_intact.subtract(orf1_orf2_intact, A=True),
    "Full length corrupted": full_length.subtract(orf1_orf2_intact, A=True).subtract(orf2_intact, A=True)
}
l1s = {k: v.subtract(blacklisted, A=True) for k, v in l1s.items()}

flag = BedTool(utils.paths.FLAG_CURAX_14H).sort().subtract(blacklisted, A=True)
coverage_flag = {l1: len(regions.intersect(flag, u=True, wa=True)) / len(regions) for l1, regions in l1s.items()}

z22 = BedTool(utils.paths.Z22_CURAX_14H).sort().subtract(blacklisted, A=True)
coverage_z22 = {l1: len(regions.intersect(z22, u=True, wa=True)) / len(regions) for l1, regions in l1s.items()}
########################################################################################################################

colors = {
    "ORF1 & ORF2 intact": "#FF7888",
    "ORF2 intact": "#FFCED6",
    "Full length corrupted": "#D5EEE1"
}
order = ["ORF1 & ORF2 intact", "ORF2 intact", "Full length corrupted"]

barw = 0.1
gap = barw
pad = barw / 2

anno_vertical_offset = 0.7
yticks = [5, 10, 15, 20, 25]
lw = 0.75
fontsize = "large"

########################################################################################################################

fig = plt.figure()
ax = fig.gca()
ax.set_xticks([])

ax.set_yticks(yticks)
ax.set_yticklabels([f"{y}%" for y in yticks], fontsize=fontsize)

# Z22
offset = pad
for ind, k in enumerate(order):
    value = round(coverage_z22[k] * 100, 1)
    rectangles = ax.bar(offset + ind * barw, value, color=colors[k], width=barw, edgecolor='black', align='edge')
    utils.miscellaneous.annotate_rectangles_with_values(rectangles, ax)
ax.text(offset + len(order) * barw / 2, -anno_vertical_offset, "Z22 peaks",
        horizontalalignment='center', verticalalignment='top', fontsize=fontsize)

# FLAG
offset = offset + len(order) * barw + gap
for ind, k in enumerate(order):
    value = round(coverage_flag[k] * 100, 1)
    rectangles = ax.bar(offset + ind * barw, value, color=colors[k], width=barw,
                        linewidth=lw, edgecolor='black', align='edge')
    utils.miscellaneous.annotate_rectangles_with_values(rectangles, ax)
ax.text(offset + len(order) * barw / 2, -anno_vertical_offset, "FLAG peaks",
        horizontalalignment='center', verticalalignment='top', fontsize=fontsize)

ax.set_ylim((0, yticks[-1] * 1.1))
ax.set_xlim(0, offset + len(order) * barw + pad)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)

legend = [mpatches.Patch(facecolor=colors[k], label=k, edgecolor='black', linewidth=lw) for k in order]
ax.legend(loc='upper right', handles=legend, frameon=False, fontsize=fontsize)
ax.set_title("LINE-1s covered by curaxin-induced peaks", fontsize=fontsize).set_position([.5, 1.05])

saveto = pathlib.Path(__file__).name.replace(".py", ".svg")
fig.savefig(utils.paths.RESULTS.joinpath(saveto), bbox_inches="tight", pad_inches=0)
