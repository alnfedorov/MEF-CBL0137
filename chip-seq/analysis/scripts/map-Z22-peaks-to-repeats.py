import pathlib
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
from HTSeq import GenomicInterval, GenomicArrayOfSets
from pybedtools import BedTool
from tqdm import tqdm

import utils
from utils.miscellaneous import is_valid_chromosome

########################################################################################################################
blacklisted = BedTool(utils.paths.BLACKLIST).sort()
repmasker = BedTool(utils.paths.REPEATMASKER).filter(is_valid_chromosome).sort().subtract(blacklisted)

# Build genomic index for repeats
genindex = GenomicArrayOfSets('auto', stranded=False)
for repeat in tqdm(repmasker):
    name, family, cls = utils.repmasker.classify(repeat.name)
    cls = cls.replace("?", "")

    interval = GenomicInterval(repeat.chrom, repeat.start, repeat.end)
    if name in ("L1Md_A", "L1Md_T"):
        genindex[interval] += name
    elif family == "L1":
        genindex[interval] += 'other LINE-1s'
    elif cls == "LINE":
        genindex[interval] += 'other LINE'
    else:
        genindex[interval] += cls

########################################################################################################################
# Map each peak to repeats
z22 = BedTool(utils.paths.Z22_CURAX_14H).filter(is_valid_chromosome).sort()
z22 = z22.subtract(blacklisted, A=True)

coverage = defaultdict(float)
for peak in tqdm(z22):
    interval = GenomicInterval(peak.chrom, peak.start, peak.end)

    # gather overlapping peaks
    peakcoverage = defaultdict(int)
    for iv, val in genindex[interval].steps():
        for v in val:
            peakcoverage[v] += iv.length

    # calculate length not covered by repeats
    totalinter = sum(peakcoverage.values())
    if totalinter > interval.length:  # rare event when repeats overlap
        peakcoverage['repeat-free'] = 0
    else:
        peakcoverage['repeat-free'] = interval.length - totalinter
        assert (interval.length - totalinter) >= 0

    # Total weight = 1, each class gets a fraction proportional to the intersection
    totalinter = sum(peakcoverage.values())
    for k, v in peakcoverage.items():
        coverage[k] += v / totalinter

# get percentages
totalpeaks = sum(coverage.values())
coverage = {k: v / totalpeaks * 100 for k, v in coverage.items()}

########################################################################################################################

# Merge rare categories
threshold = 3
coverage['other repeats'] = 0
for k in list(coverage.keys()):
    if coverage[k] < threshold and k not in ("L1Md_A", "L1Md_T", "blacklisted"):
        coverage['other repeats'] += coverage.pop(k)

# rename for the convenience
coverage['simple_repeats'] = coverage.pop('Simple_repeat')

########################################################################################################################
order = ['other LINE-1s', 'L1Md_A', 'L1Md_T', 'repeat-free', 'other repeats', 'LTR', 'simple_repeats']
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

fig = plt.figure()
ax = fig.gca()

colors = [colorscheme[k] for k in repeats]
wedges, _ = ax.pie(shares, startangle=0, wedgeprops=dict(width=0.5), colors=colors)

kw = dict(arrowprops=dict(arrowstyle="-", relpos=(0, 0.5)), va="center")

for i, p in enumerate(wedges):
    p.set_linewidth(0.3)
    p.set_edgecolor('white')

    ang = (p.theta2 - p.theta1) / 2. + p.theta1
    y = np.sin(np.deg2rad(ang))
    x = np.cos(np.deg2rad(ang))

    xytext = (x * 1.2, y * 1.2)
    xy = (x * 0.92, y * 0.92)

    horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
    ax.annotate(repeats[i], xy=xy, xytext=xytext,
                horizontalalignment=horizontalalignment, **kw, fontsize=fontsize)

title = ax.set_title("ChIP-seq by Z22 Ab\n(CBL0137-14h)", horizontalalignment='center', fontsize=fontsize)

saveto = pathlib.Path(__file__).name.replace(".py", ".svg")
fig.savefig(utils.paths.RESULTS.joinpath(saveto), bbox_inches="tight", pad_inches=0)
