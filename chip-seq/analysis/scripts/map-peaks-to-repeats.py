import pathlib
from collections import defaultdict
from typing import Dict

import matplotlib.pyplot as plt
from HTSeq import GenomicInterval, GenomicArrayOfSets
from matplotlib import patches as mpatches
from pybedtools import BedTool
from tqdm import tqdm

import utils
from utils.miscellaneous import is_valid_chromosome


def mappeaks(bed: BedTool, index: GenomicArrayOfSets) -> Dict[str, float]:
    coverage = defaultdict(float)
    for peak in tqdm(bed):
        interval = GenomicInterval(peak.chrom, peak.start, peak.end)

        # gather overlapping peaks
        peakcoverage = defaultdict(int)
        for iv, val in index[interval].steps():
            for v in val:
                peakcoverage[v] += iv.length

        # calculate length not covered by repeats
        totalinter = sum(peakcoverage.values())
        if totalinter > interval.length:  # rare event when repeats overlap
            peakcoverage['Repeat-free'] = 0
        else:
            peakcoverage['Repeat-free'] = interval.length - totalinter
            assert (interval.length - totalinter) >= 0

        # Total weight = 1, each class gets a fraction proportional to the intersection
        totalinter = sum(peakcoverage.values())
        for k, v in peakcoverage.items():
            coverage[k] += v / totalinter

    # get percentages
    totalpeaks = sum(coverage.values())
    coverage = {k: v / totalpeaks * 100 for k, v in coverage.items()}

    # Merge rare categories
    threshold = 3
    coverage['Other repeats'] = 0
    for k in list(coverage.keys()):
        if coverage[k] < threshold and k not in ("L1Md_A", "L1Md_T", "Simple_repeat", "Blacklisted"):
            coverage['Other repeats'] += coverage.pop(k)
    return coverage


def makeplot(coverage: Dict[str, float], antibody: str):
    # rename for the convenience
    coverage['Simple_repeats'] = coverage.pop('Simple_repeat', 0)

    order = ['Other LINE-1s', 'L1Md_A', 'L1Md_T', 'Repeat-free', 'Other repeats', 'LTR', 'Simple_repeats']
    fontsize = 'large'

    repeats = [k.replace("_", " ") for k in order]
    shares = [coverage[k] for k in order]

    colorscheme = {
        "Blacklisted": "#717BF9",
        "Repeat-free": "#C3E9FE",
        "Other LINE-1s": "#FEF889",
        "LTR": "#7B94E6",
        "SINE": "#186E9D",
        "Simple repeats": "#BCE98B",
        "Other repeats": "#75C66B",

        "L1Md A": "#FE9D5D",
        "L1Md T": "#EF4B57"
    }

    fig = plt.figure()
    ax = fig.gca()

    colors = [colorscheme[k] for k in repeats]
    wedges, _ = ax.pie(shares, startangle=0, wedgeprops=dict(width=0.5, edgecolor='#585A5B'), colors=colors)

    ax.set_title(f"ChIP-seq by {antibody} Ab\n(CBL0137-14h)", horizontalalignment='center', fontsize=fontsize)

    legend = [mpatches.Patch(facecolor=colorscheme[k], label=k, edgecolor='black', linewidth=0.75) for k in repeats]
    ax.set_xlim(-1, 3)
    ax.legend(handles=legend, frameon=False, fontsize='x-large', loc='right')

    saveto = pathlib.Path(__file__).name.replace(".py", f".{antibody}.eps")
    fig.savefig(utils.paths.RESULTS.joinpath(saveto), bbox_inches="tight", pad_inches=0)


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
        genindex[interval] += 'Other LINE-1s'
    elif cls == "LINE":
        genindex[interval] += 'Other LINE'
    else:
        genindex[interval] += cls

# build plots
for name, bed in {"Z22": utils.paths.Z22_CURAX_14H, "FLAG": utils.paths.FLAG_CURAX_14H}.items():
    print(name)
    bed = BedTool(bed).filter(is_valid_chromosome).sort()
    bed = bed.subtract(blacklisted, A=True)
    coverage = mappeaks(bed, genindex)
    makeplot(coverage, name)
