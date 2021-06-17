from collections import defaultdict
from functools import lru_cache

import numpy as np
import pandas as pd
import tqdm
from HTSeq import GenomicInterval, GenomicArrayOfSets
from pybedtools import BedTool

import utils


@lru_cache(maxsize=1)
def build_genes_index() -> GenomicArrayOfSets:
    INDEX = GenomicArrayOfSets('auto', stranded=False)
    for region, bed in [("utr3", utils.paths.GENCODE_3UTR), ("utr5", utils.paths.GENCODE_5UTR),
                        ("exons", utils.paths.GENCODE_EXONS), ("introns", utils.paths.GENCODE_INTRONS)]:
        for interval in BedTool(bed):
            transcript_id = interval.name.split("_")[0]
            gname = utils.ensembl.transcript_to_gene_name(transcript_id)
            interval = GenomicInterval(interval.chrom, interval.start, interval.end)
            INDEX[interval] += (region, gname)
    return INDEX


def map_to_genes(df: pd.DataFrame):
    index = build_genes_index()

    gene = []
    geneloc = []
    for coord in tqdm.tqdm(df['Coordinates']):
        chrom, buffer = coord.split(":")
        start, end = buffer.split("-")
        coord = GenomicInterval(chrom, int(start), int(end))

        keys = set()
        for iv, intersected in index[coord].steps():
            if iv.length > 0:
                keys |= intersected
        features = defaultdict(list)
        for key in keys:
            features[key[0]].append(key[1])

        gname, gloc = None, None
        for region in ["utr3", "utr5", "introns", "exons"]:
            if region not in features:
                continue
            gname = ";".join(features[region])
            gloc = region
            break
        gene.append(gname if gname else "na")
        geneloc.append(gloc if gloc else "intergenic")
    df['Gene'] = gene
    df['Location'] = geneloc
    return df


def loadediting(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df['Coordinates'] = [f'{chrom}:{start}-{end}'
                         for chrom, start, end in zip(df['GenomicRegion'], df['Start'], df['End'])]

    df = df[df['A2G'] != "Unknown"]

    df['Mismatches'] = np.where(df['A2G'] == "+", df['NumOfA2GMismatches'], df['NumOfT2CMismatches'])
    df['Matches'] = np.where(df['A2G'] == "+", df['NumOfA'], df['NumOfT'])
    df['EditingIndex'] = df['Mismatches'] / (df['Mismatches'] + df['Matches']) * 100

    a2g = df['NumOfA2GMismatches'] / (df['NumOfA2GMismatches'] + df['NumOfA']) * 100
    t2c = df['NumOfT2CMismatches'] / (df['NumOfT2CMismatches'] + df['NumOfT']) * 100
    index = np.where(df['A2G'] == '+', a2g, t2c)
    diff = np.abs(df['EditingIndex'] - index)
    assert diff.max() < 1e-3, diff.max()

    # Keep repeats, where each A/T is covered by at least X reads
    reads_per_A = df['TotalCoverageAtAPositions'] / df['NumOfAPositionsCovered']
    reads_per_T = df['TotalCoverageAtTPositions'] / df['NumOfTPositionsCovered']
    df['Gene'] = np.where(df['A2G'] == "+", df['SenseGeneCommonName'], df['AntisenseGeneCommonName'])
    df['MeanCoveragePerSite'] = np.where(df['A2G'] == "+", reads_per_A.round(2), reads_per_T.round(2))
    df['InferredStrand'] = df['A2G']
    df = df[[
        'Coordinates', 'InferredStrand', 'Gene', 'EditingIndex', 'MeanCoveragePerSite', 'Mismatches', 'Matches'
    ]]
    return df.set_index(['Coordinates', 'InferredStrand', 'Gene'])


def loadsample(sample: str, minreads_per_site: float = 10, editingthr: float = 0.5):
    perrepdf = []
    for folder in utils.paths.EDITING_INDEX_VALUES.iterdir():
        file = folder.joinpath(sample, "StrandDerivingCountsPerRegion.csv")
        df = loadediting(file).reset_index()
        df = df[(df['MeanCoveragePerSite'] > minreads_per_site) & (df['EditingIndex'] > editingthr)]
        df = map_to_genes(df)
        df['RepeatCls'] = folder.name
        perrepdf.append(df)
    df = pd.concat(perrepdf)
    return df


df = loadsample("RIP-ADAR1-WT-IFNb-72h-Z22")

ISG = utils.ISG.names()
df['ISG'] = df['Gene'].apply(lambda genes: any(x in ISG for x in genes.split(";")) if genes != "na" else False)

df['Count'] = 0
df = df.groupby(['ISG', 'RepeatCls', 'Location'])['Count'].count().reset_index()
df = df.pivot(['ISG', 'RepeatCls'], 'Location', 'Count').fillna(0).astype(int)
df['GrandTotal'] = df['exons'] + df['intergenic'] + df['introns'] + df['utr3'] + df['utr5']

df.to_csv(
    utils.paths.RESULTS.joinpath("Editing-distribution(Z22 & editing index > 0.5% & mean coverage > 10).csv")
)
