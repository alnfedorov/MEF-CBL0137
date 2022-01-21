import pathlib

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd

import utils

plt.rcParams['svg.fonttype'] = 'none'

df = pd.read_csv(utils.paths.EDITING_INDEX_VALUES.joinpath("SINE", "EditingIndex.csv"))
df = df.set_index('Sample')['A2GEditingIndex']

shares = {
    "ADAR1-WT": {
        "IFNβ (-)": df["RIP-ADAR1-WT-IFNb-0h-Z22"],
        "IFNβ (+)": df["RIP-ADAR1-WT-IFNb-72h-Z22"]
    },
    # "ADAR1-KO": {
    #     "IFNβ (-)": df["RIP-ADAR1-KO-IFNb-0h-Z22"],
    #     "IFNβ (+)": df["RIP-ADAR1-KO-IFNb-72h-Z22"]
    # },
    "ADAR1/ZBP1 dKO": {
        "IFNβ (-)": df["ADAR1-KO+ZBP1-KO_IFNb-0h_Z22"],
        "IFNβ (+)": df["ADAR1-KO+ZBP1-KO_IFNb-48h_Z22"]
    }
}

colors = {"ADAR1-WT": "#EF553B", "ADAR1/ZBP1 dKO": "#636EFA"}  # , "ADAR1-KO": "#666666"}

order = ["IFNβ (-)", "IFNβ (+)"]
elemorder = ["ADAR1/ZBP1 dKO", "ADAR1-WT"]
# elemorder = ["ADAR1/ZBP1 dKO", "ADAR1-KO", "ADAR1-WT"]

barw = 0.5
gap = barw / 4
pad = barw / 4

anno_vertical_offset = 0.02
yticks = [0.2, 0.4, 0.6, 0.8, 1.0]
lw = 0.75
fontsize = "large"
fontsizelbl = "x-large"

fig = plt.figure(figsize=(7, 6))
ax = fig.gca()
ax.set_xticks([])

ax.set_yticks(yticks)
ax.set_ylabel('Alu editing index (%)', fontsize=fontsizelbl)
ax.set_yticklabels([f"{y}" for y in yticks], fontsize=fontsize)

offset = pad
for antibody in order:
    for ind, celltype in enumerate(elemorder):
        value = round(shares[celltype][antibody], 3)
        rectangles = ax.bar(offset + ind * barw, value, color=colors[celltype],
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
ax.set_title(f"RNAs bound by Z22 Ab", fontsize=fontsizelbl).set_position([.5, 1.05])

saveto = pathlib.Path(__file__).name.replace(".py", f".svg")
fig.savefig(utils.paths.RESULTS.joinpath(saveto), bbox_inches="tight", pad_inches=0)

saveto = pathlib.Path(__file__).name.replace(".py", f".eps")
fig.savefig(utils.paths.RESULTS.joinpath(saveto), bbox_inches="tight", pad_inches=0)
