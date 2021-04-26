import pathlib
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pybedtools
from pybedtools import BedTool
from tqdm import tqdm

import utils
from utils.miscellaneous import is_valid_chromosome

########################################################################################################################
blacklisted = BedTool(utils.paths.BLACKLIST).filter(is_valid_chromosome).sort()
repmasker = BedTool(utils.paths.REPEATMASKER).filter(is_valid_chromosome).sort().subtract(blacklisted)

# calculate coverage for each repeat family
coverage = defaultdict(int)
coverage['blacklisted'] = blacklisted.total_coverage()
for repeat in tqdm(repmasker):
    name, family, cls = utils.repmasker.classify(repeat.name)
    cls = cls.replace("?", "")
    if name in ("L1Md_A", "L1Md_T"):
        coverage[name] += repeat.length
    elif family == "L1":
        coverage['other LINE-1s'] += repeat.length
    elif cls == "LINE":
        coverage['other LINE'] += repeat.length
    else:
        coverage[cls] += repeat.length

# calculate repeat-free fraction of the genome
mm10 = {k: v for k, (_, v) in pybedtools.chromsizes('mm10').items() if is_valid_chromosome(k)}
mm10_totalsize = sum(mm10.values())

coverage['repeat-free'] = mm10_totalsize - sum(coverage.values())
assert sum(coverage.values()) == mm10_totalsize

# normalize coverage and get percents
coverage = {k: v / mm10_totalsize * 100 for k, v in coverage.items()}
coverage['other repeats'] = coverage.pop('Other')

# merge rare repeats
threshold = 2
for k in list(coverage.keys()):
    if coverage[k] < threshold and k not in ("L1Md_A", "L1Md_T", "blacklisted"):
        coverage['other repeats'] += coverage.pop(k)

# rename for simplicity
coverage['simple_repeats'] = coverage.pop('Simple_repeat')

########################################################################################################################
# make the plot
order = ['repeat-free', 'blacklisted', 'other repeats', 'SINE', 'LTR',
         'simple_repeats', 'other LINE-1s', 'L1Md_A', 'L1Md_T']
fontsize = 'large'

repeats = [k.replace("_", " ") for k in order]
shares = [coverage[k] for k in order]

colorscheme = {
    "blacklisted": "#7f7f7f",
    "repeat-free": "#2CA674",
    "other LINE-1s": "#F49E2D",
    "LTR": "#D079B1",
    "SINE": "#882255",
    "simple repeats": "#4B0092",
    "other repeats": "#636EFA",

    "L1Md A": "#E75C21",
    "L1Md T": "#F6E358"
}

fig = plt.figure(figsize=(16, 8))
ax = fig.gca()

colors = [colorscheme[k] for k in repeats]
wedges, _ = ax.pie(shares, startangle=90, wedgeprops=dict(width=0.5), colors=colors)

kw = dict(arrowprops=dict(arrowstyle="-", relpos=(0, 0.5)), va="center")

for i, p in enumerate(wedges):
    p.set_linewidth(0.3)
    p.set_edgecolor('white')
    if repeats[i] in ("L1Md A", "L1Md T"):
        continue

    ang = (p.theta2 - p.theta1) / 2. + p.theta1
    y = np.sin(np.deg2rad(ang))
    x = np.cos(np.deg2rad(ang))

    xytext = (x * 1.2, y * 1.2)
    xy = (x * 0.92, y * 0.92)

    horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
    ax.annotate(repeats[i], xy=xy, xytext=xytext,
                horizontalalignment=horizontalalignment, **kw, fontsize=fontsize)

kw['va'] = 'bottom'

# L1Md A
i = repeats.index("L1Md A")

ang = (wedges[i].theta2 - wedges[i].theta1) / 2. + wedges[i].theta1
y = np.sin(np.deg2rad(ang))
x = np.cos(np.deg2rad(ang))

xytext = (x * 1.08, y * 1.08)
xy = (x * 0.92, y * 0.92)
ax.annotate(repeats[i], xy=xy, xytext=xytext,
            horizontalalignment='left', **kw, fontsize=fontsize)

# L1Md T
i = repeats.index("L1Md T")

ang = (wedges[i].theta2 - wedges[i].theta1) / 2. + wedges[i].theta1
y = np.sin(np.deg2rad(ang))
x = np.cos(np.deg2rad(ang))

xytext = (x * 1.25, y * 1.25)
xy = (x * 0.92, y * 0.92)

ax.annotate(repeats[i], xy=xy, xytext=xytext, horizontalalignment='left', **kw, fontsize=fontsize)

ax.set_title("Mouse genome coverage by repeats", fontsize=fontsize).set_position([.5, 1.05])

saveto = pathlib.Path(__file__).name.replace(".py", ".svg")
fig.savefig(utils.paths.RESULTS.joinpath(saveto), bbox_inches="tight", pad_inches=0)
