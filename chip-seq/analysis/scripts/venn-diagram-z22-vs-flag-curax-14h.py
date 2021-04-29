import pathlib

import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib_venn._venn2 import *
from pybedtools import BedTool

import utils


def intersections(peaks):
    peaks = {k: v.sort() for k, v in peaks.items()}
    subsets = {}
    subsets["11"] = peaks["10"].intersect(peaks["01"], u=True, wa=True)

    subsets["10"] = peaks["10"].subtract(subsets["11"], A=True)
    subsets["01"] = peaks["01"].subtract(subsets["11"], A=True)

    subsets = {k: len(v) for k, v in subsets.items()}
    if subsets["11"] + subsets["10"] != len(peaks["10"]):
        subsets["10"] += len(peaks["10"]) - (subsets["11"] + subsets["10"])
    if subsets["11"] + subsets["01"] != len(peaks["01"]):
        subsets["01"] += len(peaks["01"]) - (subsets["11"] + subsets["01"])

    assert subsets["11"] + subsets["10"] == len(peaks["10"])
    assert subsets["11"] + subsets["01"] == len(peaks["01"])

    return subsets


blacklisted = BedTool(utils.paths.BLACKLIST)
z22 = BedTool(utils.paths.Z22_CURAX_14H).subtract(blacklisted, A=True)
flag = BedTool(utils.paths.FLAG_CURAX_14H).subtract(blacklisted, A=True)

l1s = BedTool(utils.paths.REPEATMASKER).filter(
    lambda x: utils.repmasker.classify(x.name)[1] == "L1"
)

# Calculate total intact l1s inside the intersection
intactl1s = BedTool(utils.paths.L1BASE_ORF1_ORF2_INTACT).sort()
intactl1s = l1s.intersect(intactl1s, u=True, wa=True)
intersectedl1s = len(z22.intersect(flag, u=True, wa=True).intersect(intactl1s, u=True, wa=True))

subsets = intersections({"01": z22, "10": flag})

########################################################################################################################
# make the plot
colors = ("#7CD7F5", "#FF4474", "#7CC0E4")
alphas = (1, 1, 1)
fontsize = 'x-large'

fig = plt.figure()
ax = fig.gca()

areas = compute_venn2_areas([subsets.get(t, 0) for t in ['10', '01', '11']])
centers, radii = solve_venn2_circles(areas)

regions = compute_venn2_regions(centers, radii)
prepare_venn_axes(ax, centers, radii)

# Create and add patches and subset labels
for r, c, alpha in zip(regions, colors, alphas):
    p = r.make_patch()
    p.set_facecolor(c)
    p.set_alpha(alpha)
    # p.set_edgecolor((0, 0, 0, 0))
    ax.add_patch(p)

for c, r in zip(centers, radii):
    ax.add_patch(patches.Circle(c, r, linewidth=1.5, fill=False, color='k'))

kw = dict(arrowprops=dict(arrowstyle="-", relpos=(0.5, 1)), va="center", ha="center")
ax.annotate("FLAG peaks", (-0.46, -0.22), (-0.57, -0.42), **kw, fontsize=fontsize)

kw['arrowprops']['relpos'] = (0, 0.5)
kw['ha'] = 'left'
ax.annotate("Z22 peaks", (0.58, 0), (0.7, 0), **kw, fontsize='x-large')

z22_unique, flag_unique, intersection = subsets['01'], subsets['10'], subsets['11']

ax.text(-0.12, 0, f"{intersection}", fontstyle='oblique', va='center', ha='center', fontsize=fontsize)
ax.text(0.4, 0, f"{z22_unique}", fontstyle='oblique', va='center', ha='center', fontsize=fontsize)
ax.text(-0.49, 0, f"{flag_unique}", fontstyle='oblique', va='center', ha='center', fontsize=fontsize)

kw = dict(arrowprops=dict(arrowstyle="->", relpos=(0.5, 1)), va="center", ha="center")
ax.annotate(f"{intersectedl1s}({round(intersectedl1s / intersection * 100)}%) "
            f"shared peaks map to ORF1 & ORF2 intact LINE-1s",
            (-0.12, -0.03), (-0.0145, -0.607), fontsize='large', **kw)

ax.set_title("Z22 ChIP-seq (curaxin 14h)\nFLAG-ZBP1(PLUS)", horizontalalignment='center', fontsize=fontsize)
fig.tight_layout()
saveto = pathlib.Path(__file__).name.replace(".py", ".svg")
fig.savefig(utils.paths.RESULTS.joinpath(saveto), bbox_inches="tight", pad_inches=0)
