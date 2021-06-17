import pickle
from collections import defaultdict

import matplotlib.pyplot as plt
import pandas as pd

import utils


def infermeta(sample: str) -> str:
    if "IgG" in sample:
        antibody = "IgG"
    elif "Z22" in sample:
        antibody = "Z22"
    elif "J2" in sample:
        antibody = "J2"
    else:
        antibody = pd.NA

    if "0h" in sample:
        treatment = "IFNb-0h"
    else:
        assert "72h" in sample
        treatment = "IFNb-72h"

    if "KO" in sample:
        cells = "KO"
    else:
        assert "WT" in sample
        cells = "WT"
    return antibody, treatment, cells


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
        mappings[key] = "ambiguous"
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
df['antibody'], df['treatment'], df['cells'] = zip(*df['index'].apply(infermeta))
df = df[~df['antibody'].isna()].drop(columns=['index'])
df = df.rename(columns={
    "protein_coding": "mRNA", "other_RNA": "Other RNA", "ambiguous": "Ambiguous", "intron": "Intronic",
    "intergenic": "Intergenic"
})

# make plots
colorsheme = {
    'mRNA': '#FB8A3A', 'miRNA': '#EA555A', 'Other RNA': '#509F55', 'Ambiguous': '#B2799F',
    'Intronic': '#447AA4', 'Intergenic': '#666666'
}

# stacked barplot
fig, axes = plt.subplots(2, 2, sharex=True, figsize=(16, 16))
axes = axes.ravel()

df = df.sort_values(by=['antibody']).set_index('antibody')
locations = ['mRNA', 'miRNA', 'Other RNA', 'Ambiguous', 'Intronic', 'Intergenic']
for treatment in ['IFNb-0h', 'IFNb-72h']:
    for cells in ['WT', 'KO']:
        axes[0].set_title(f"{cells} - {treatment}")
        df[(df['cells'] == cells) & (df['treatment'] == treatment)][locations].plot.barh(
            stacked=True, ax=axes[0], color=colorsheme
        )
        axes = axes[1:]
df = df.reset_index()

fig.savefig(utils.paths.RESULTS.joinpath("fragments-per-genomic-features.svg"), bbox_inches="tight", pad_inches=0)

# pie-charts

