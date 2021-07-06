import numpy as np
import pandas as pd
import scipy.stats as stats
from HTSeq import GenomicArrayOfSets, GenomicInterval
from pybedtools import BedTool
from tqdm import tqdm

import utils


# fisher exact utils
def expected_counts(table: np.ndarray) -> np.ndarray:
    coltotal, rowstotal, grandtotal = table.sum(axis=0), table.sum(axis=1), table.sum()
    expected = np.asarray([
        [rowstotal[0] * coltotal[0], rowstotal[0] * coltotal[1]],
        [rowstotal[1] * coltotal[0], rowstotal[1] * coltotal[1]]
    ]) / grandtotal
    return expected


def print_table(table: np.ndarray, rows, columns):
    print(f"{'':<20}{columns[0]:^20}{columns[1]:^20}")
    for row, name in zip(table, rows):
        print(f"{name:<20}", end="")
        print(f"{row[0]:^20}{row[1]:^20}")


# Only a subset (%?) of Z22-enriched mRNAs had SINEs in their 3’UTRs.

# gather Z22 enriched genes
Z22enriched = {}
for cell in ["ADAR1-WT", "ADAR1-KO"]:
    for ifnb in ["0h", "72h"]:
        df = pd.read_csv(utils.paths.DE_GENES.joinpath(f"{cell}-{ifnb}-Z22-vs-IgG.csv"))
        genes = set(df[(df['padj'] < 0.1) & (df['log2FoldChange'] > 0.585)]['Gene name'])
        Z22enriched[f"{cell}-{ifnb}"] = genes
Z22enriched["All"] = set.union(*Z22enriched.values())

# gather labeled 3`UTRs
index = GenomicArrayOfSets('auto', stranded=False)
for interval in tqdm(BedTool(utils.paths.GENCODE_3UTR)):
    transcript_id = interval.name.split("_")[0]
    gname = utils.ensembl.transcript_to_gene_name(transcript_id)
    interval = GenomicInterval(interval.chrom, interval.start, interval.end)
    index[interval] += gname

# map SINEs to 3`UTRs
genes_with_sines = set()
for repeat in tqdm(BedTool(utils.paths.REPMASKER)):
    _, _, reptype = utils.repmasker.classify(repeat.name)
    if reptype != "SINE":
        continue
    interval = GenomicInterval(repeat.chrom, repeat.start, repeat.end)
    for iv, intersected in index[interval].steps():
        if iv.length > 0:
            genes_with_sines |= intersected

# calculate % of Z22-enriched genes with SINEs in 3`UTRs
print("% of Z22 enriched genes with SINEs in 3`UTRs")
for key, genes in Z22enriched.items():
    intersection = genes.intersection(genes_with_sines)
    print(f"{key:<15} => {len(intersection)}/{len(genes)}({len(intersection) / len(genes):.2f})")

# Most of the ADAR1-edited RNAs pulled down by the Z22 antibody were SINEs and simple repeats, and a notable fraction
# of these were found in ISG-encoded mRNA transcripts (Fig. S1dS1e).

# build the background (genes with repeats inside)
repeats = {"SINE": [], "Simple_repeat": []}
for interval in BedTool(utils.paths.REPMASKER):
    _, _, cls = utils.repmasker.classify(interval.name)
    if cls in repeats:
        repeats[cls].append(interval)

genes = BedTool(utils.paths.GENCODE_GENES).sort()
bckggenes = {k: genes.intersect(BedTool(v).sort(), wa=True, u=True) for k, v in repeats.items()}
bckggenes = {k: set(utils.ensembl.transcript_to_gene_name(x.name.split("_")[0]) for x in v) for k, v in
             bckggenes.items()}

ISG = utils.ISG.names()
bckgisg = {k: v.intersection(ISG) for k, v in bckggenes.items()}

df = pd.read_csv(utils.paths.RESULTS.joinpath(
    "Editing-distribution(WT & IFNb-72h & Z22 & editing index > 0.5% & mean coverage > 10).csv"
))[["ISG", "GrandTotal", "RepeatCls"]]

# These p-values are strange, R-calculated are used in the paper
for repeat in ["SINE", "Simple_repeat"]:
    allgenes = len(bckggenes[repeat])
    allisg = len(bckgisg[repeat])
    pooledISG = int(df['GrandTotal'][(df['RepeatCls'] == repeat) & (df["ISG"])])
    poolednonISG = int(df['GrandTotal'][(df['RepeatCls'] == repeat) & (~df["ISG"])])
    pooledall = pooledISG + poolednonISG
    print(repeat)
    print(
        f"(Sample) ISG among genes with edited {repeat}: {pooledISG} / {pooledall}({pooledISG / pooledall * 100:.2f}%)")
    print(f"(Background) ISG among genes with {repeat}: {allisg}/{allgenes}({allisg / allgenes * 100:.2f}%)")
    print(f"R code: phyper({pooledISG}-1, {allisg}, {allgenes - allisg}, {pooledall}, lower.tail=FALSE)")
    print(
        f"p-value: {stats.hypergeom(allgenes, allisg, poolednonISG + pooledISG).sf(pooledISG - 1)}"
    )
    print()

mask = (df['RepeatCls'] == "SINE") | (df['RepeatCls'] == "Simple_repeat")

repofinterest = df[mask].groupby("ISG")['GrandTotal'].sum().to_dict()
bckgrep = df[~mask].groupby("ISG")['GrandTotal'].sum().to_dict()

table = [
    [repofinterest[True], repofinterest[False]],
    [bckgrep[True], bckgrep[False]]
]
table = np.asarray(table)
pvalue = stats.fisher_exact(table)[1]

expected = expected_counts(table).round(2)

print(f"Fisher exact test p-value: {pvalue}")
for t in table, expected:
    print_table(t, ["SINEs/Simple", "Other"], ["ISG", "non-ISG"])
    print()

# 3’UTRs of ISG mRNAs were disproportionately targeted for editing by ADAR1 p150
df = pd.read_csv(utils.paths.RESULTS.joinpath(
    "RAW-Editing-distribution(WT & IFNb-72h & Z22 & editing index > 0.5% & mean coverage > 10).csv"
))

# 1 - map edited repeats to 3`UTRs of genes
genes_with_annotated_utr3 = set()
index = GenomicArrayOfSets('auto', stranded=False)
for interval in tqdm(BedTool(utils.paths.GENCODE_3UTR)):
    transcript_id = interval.name.split("_")[0]
    gname = utils.ensembl.transcript_to_gene_name(transcript_id)
    interval = GenomicInterval(interval.chrom, interval.start, interval.end)
    index[interval] += gname
    genes_with_annotated_utr3.add(gname)
ISG_with_annotated_utr3 = genes_with_annotated_utr3.intersection(utils.ISG.names())

genes_with_edited_utr3 = set()
for interval in df['Coordinates']:
    chrom, buffer = interval.split(":")
    start, end = buffer.split("-")
    interval = GenomicInterval(chrom, int(start), int(end))
    for iv, intersected in index[interval].steps():
        if iv.length > 0:
            genes_with_edited_utr3 |= intersected
ISG_with_edited_utr3 = genes_with_edited_utr3.intersection(utils.ISG.names())

genes_with_annotated_utr3 = len(genes_with_annotated_utr3)
ISG_with_annotated_utr3 = len(ISG_with_annotated_utr3)
genes_with_edited_utr3 = len(genes_with_edited_utr3)
ISG_with_edited_utr3 = len(ISG_with_edited_utr3)

print("Editing in 3`UTRs of ISG")
print(f"(Sample)ISG having edited 3`UTRs: {ISG_with_edited_utr3}: {ISG_with_edited_utr3} / {genes_with_edited_utr3}")
print(
    f"(Background)ISG having annotated 3`UTRs: {ISG_with_annotated_utr3}: {ISG_with_annotated_utr3} / {genes_with_annotated_utr3}")
print(
    f"p-value: {stats.hypergeom(genes_with_annotated_utr3, ISG_with_annotated_utr3, genes_with_edited_utr3).sf(ISG_with_edited_utr3 - 1)}"
)
