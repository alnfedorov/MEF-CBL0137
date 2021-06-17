import pathlib

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd

import utils

df = pd.read_csv(utils.paths.EDITING_INDEX_VALUES.joinpath("SINE", "EditingIndex.csv"))
shares = {"IgG": {}, "J2": {}, "Z22": {}}
for _, row in df.iterrows():
    if "WT" not in row['Sample']:
        continue

    meta = row['Sample'].split('-')
    antibody, ifnb = meta[-1], meta[-2]
    shares[antibody][f"IFNβ {ifnb}"] = row['A2GEditingIndex']

colors = {"IFNβ 0h": "#5B84B0", "IFNβ 72h": "#FC766A"}

elemorder = ["IFNβ 0h", "IFNβ 72h"]
order = ["IgG", "J2", "Z22"]

barw = 0.15
gap = barw / 2
pad = barw / 2

anno_vertical_offset = 0.02
yticks = [0.15, 0.30, 0.45, 0.6, 0.75, 0.9]
lw = 0.75
fontsize = "large"
fontsizelbl = "x-large"

fig = plt.figure(figsize=(7, 6))
ax = fig.gca()
ax.set_xticks([])

ax.set_yticks(yticks)
ax.set_ylabel('Alu editing index', fontsize=fontsize)
ax.set_yticklabels([f"{y}%" for y in yticks], fontsize=fontsize)

offset = pad
for antibody in order:
    for ind, condition in enumerate(elemorder):
        value = round(shares[antibody][condition], 3)
        rectangles = ax.bar(offset + ind * barw, value, color=colors[condition],
                            width=barw, edgecolor='black', align='edge')
        utils.miscellaneous.annotate_rectangles_with_values(rectangles, ax)
    ax.text(offset + len(colors) * barw / 2, -anno_vertical_offset, antibody, fontweight="bold",
            horizontalalignment='center', verticalalignment='top', fontsize=fontsizelbl)
    offset = offset + len(colors) * barw + gap

ax.set_ylim((0, yticks[-1] * 1.1))
ax.set_xlim(0, offset)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)

legend = [mpatches.Patch(facecolor=colors[k], label=k.replace("_", " "),
                         edgecolor='black', linewidth=lw) for k in elemorder]
ax.legend(loc='upper left', handles=legend, frameon=False, fontsize=fontsizelbl)
ax.set_title("AEI, ADAR1-WT", fontsize=fontsizelbl).set_position([.5, 1.05])

saveto = pathlib.Path(__file__).name.replace(".py", ".eps")
fig.savefig(utils.paths.RESULTS.joinpath(saveto), bbox_inches="tight", pad_inches=0)
