import pickle
from collections import defaultdict
from typing import Tuple

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import patches as mpatches

import utils


def infermeta(sample: str) -> Tuple[str, str, str, str]:
    cells, treatment, antibody, rep = sample.split(".")[0].split("_")
    return cells, treatment, antibody, rep


def simplify_keys(data):
    synonyms = {
        "lncRNA": ["3prime_overlapping_ncRNA", "antisense", "bidirectional_promoter_lncRNA", "lincRNA", "macro_lncRNA",
                   "non_coding", "processed_transcript", "sense_intronic", "sense_overlapping"],
        "protein_coding": ["IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_LV_gene", "IG_V_gene", "TR_C_gene",
                           "TR_J_gene", "TR_V_gene", "TR_D_gene"]
    }
    synonyms = {subv: k for k, v in synonyms.items() for subv in v}

    # simplify keys
    newdata = defaultdict(int)
    for k, v in data.items():
        if isinstance(k, str):
            k = (k,)
        newk = []
        for subk in k:
            if subk in synonyms:
                subk = synonyms[subk]
            if "pseudogene" in subk:
                subk = "pseudogene"
            newk.append(subk)
        newdata[tuple(set(newk))] += v
    data = newdata

    # filter rRNA
    data = {k: v for k, v in data.items() if not any("rRNA" in subk for subk in k)}

    # remove TEC = to be experimentally verified
    # if some key is inside intron -> remove the intron
    # remove misc_RNA if we have other alternatives
    blacklist = ["misc_RNA", "TEC", "intron"]
    newdata = defaultdict(int)
    for k, v in data.items():
        for bl in blacklist:
            if len(k) < 2:
                break
            k = tuple(subk for subk in k if subk != bl)

        newdata[k] += v
    data = newdata

    # hard fix for ribozymes
    for k in ('ribozyme', 'lncRNA'), ('lncRNA', 'ribozyme'):
        if k in data:
            data[('ribozyme',)] += data.pop(k)
    return data


total = defaultdict(int)
samples = {}
for file in utils.paths.FRAGMENTS_TO_FEATURES.iterdir():
    with open(file.absolute(), 'rb') as stream:
        counts = pickle.load(stream)
        counts = simplify_keys(counts)
        for k, v in counts.items():
            total[k] += v
        samples[file.name] = counts

allfragments = sum(total.values())
threshold = 0.01
mappings = {}
for key, counts in total.items():
    if "protein_coding" in key and len(key) >= 2:
        mappings[key] = "other_RNA"
        continue

    if counts / allfragments > threshold:
        assert len(key) == 1
        mappings[key] = key[0]
        continue

    mappings[key] = "other_RNA"

# hard fix for misc_RNA
mappings['misc_RNA'] = "other_RNA"

newsamples = {}
for sample, counts in samples.items():
    mergedcounts = defaultdict(int)
    for k, v in counts.items():
        mergedcounts[mappings[k]] += v

    allfragments = sum(mergedcounts.values())
    newsamples[sample] = {k: v / allfragments * 100 for k, v in mergedcounts.items()}

df = pd.DataFrame(newsamples).T.reset_index()
df["cells"], df["treatment"], df["antibody"], df["rep"] = zip(*df['index'].apply(infermeta))
df = df.rename(columns={
    "protein_coding": "mRNA", "other_RNA": "Other RNA", "intron": "Intronic",
    "intergenic": "Intergenic"
})

# make plots
locations = ['mRNA', 'Intronic', 'Intergenic', 'Other RNA']
colorsheme = {
    'mRNA': '#FA954D', 'Other RNA': '#60A865',
    'Intronic': '#5586AC', 'Intergenic': '#747474'
}

# pie-charts
colors = [colorsheme[loc] for loc in locations]
for ifnb in ["IFNb-0h", "IFNb-24h", "IFNb-48h"]:
    shares = df[
        (df['treatment'] == ifnb) &
        (df['antibody'] == "Z22") &
        (df['cells'] == "ADAR1-KO+ZBP1-KO")
        ][locations]
    shares = shares.median(axis=0).to_dict()
    print(ifnb, shares['mRNA'])
    shares = [shares[loc] for loc in locations]

    fig = plt.figure()
    ax = plt.gca()
    ax.set_title(f"RIP-seq by Z22 Ab [{ifnb.replace('b', 'β')}]\n(ADAR1-KO & ZBP1-KO)",
                 horizontalalignment='center', fontsize="xx-large")

    wedges, _ = ax.pie(shares, startangle=90, colors=colors,
                       wedgeprops=dict(width=0.5, linewidth=0.5, edgecolor='#585A5B'))
    ax.set_xlim(-1, 3)
    legend = [mpatches.Patch(facecolor=colorsheme[k], label=k, edgecolor='black', linewidth=0.75) for k in locations]
    ax.legend(handles=legend, frameon=False, fontsize='x-large', loc='right')

    saveto = utils.paths.RESULTS.joinpath(f"fragments-to-genomic-features[{ifnb}].eps")
    fig.savefig(saveto, bbox_inches="tight", pad_inches=0)
