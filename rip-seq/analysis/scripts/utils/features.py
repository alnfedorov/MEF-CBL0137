from collections import defaultdict
from functools import lru_cache
from itertools import chain
from typing import Set, Dict, List, Optional, Tuple

import pandas as pd
from Bio.Seq import Seq
from HTSeq import GenomicArrayOfSets, GenomicInterval
from pybedtools import BedTool, Interval

import utils.DEG
from . import paths
from .miscellaneous import hashinter, distance


def _all_repeats() -> List[Interval]:
    return list(BedTool(paths.REPMASKER))


@lru_cache(maxsize=1)
def _repeats_index() -> GenomicArrayOfSets:
    INDEX = GenomicArrayOfSets('auto', stranded=False)
    for interval in _all_repeats():
        assert interval.strand != ".", f"{interval.name} {interval.chrom}:{interval.start}-{interval.end}"
        index = GenomicInterval(interval.chrom, interval.start, interval.end)
        INDEX[index] += interval
    return INDEX


def _load_intervals(transcripts: Set[str], files: List[str]) -> Dict[str, List[Interval]]:
    saveto = [defaultdict(list) for _ in files]
    for bed, s in zip(files, saveto):
        for record in BedTool(bed):
            tid = record.name.split("_")[0].split(".")[0]
            if tid not in transcripts:
                continue
            record.name = tid
            s[tid].append(record)

    # Assert consistency between files
    annotated = [set(x.keys()) for x in saveto]
    for i in range(len(saveto)):
        for j in range(i + 1, len(saveto)):
            for tid in annotated[i] & annotated[j]:
                irec, jrec = saveto[i][tid], saveto[j][tid]
                assert len(irec) == len(jrec)
                irec, jrec = set((x.chrom, x.start, x.end, x.strand) for x in irec), \
                             set((x.chrom, x.start, x.end, x.strand) for x in jrec)
                assert irec == jrec

    result = saveto[0]
    for records in saveto[1:]:
        for tid, items in records.items():
            if tid not in result:
                result[tid] = items
    return result


def utr3(transcripts: Set[str]) -> Dict[str, List[Interval]]:
    return _load_intervals(transcripts, [utils.paths.GENCODE_3UTR,
                                         utils.paths.RESOURCES.joinpath("Mus_musculus.GRCm38.102.UTR3.bed.gz")])


def cds(transcripts: Set[str]) -> Dict[str, List[Interval]]:
    return _load_intervals(transcripts, [utils.paths.RESOURCES.joinpath("Mus_musculus.GRCm38.102.CDS.bed.gz")])


def utr5(transcripts: Set[str]) -> Dict[str, List[Interval]]:
    return _load_intervals(transcripts, [utils.paths.GENCODE_5UTR,
                                         utils.paths.RESOURCES.joinpath("Mus_musculus.GRCm38.102.UTR5.bed.gz")])


def exons(transcripts: Set[str]) -> Dict[str, List[Interval]]:
    return _load_intervals(transcripts, [utils.paths.GENCODE_EXONS])


def mRNA(transcripts: Set[str]) -> Dict[str, List[Interval]]:
    features = {
        "3`UTR": utils.features.utr3(transcripts),
        "5`UTR": utils.features.utr5(transcripts),
        "CDS": utils.features.cds(transcripts)
    }
    exons = defaultdict(list)
    for exon in features.values():
        for tid, r in exon.items():
            exons[tid].extend(r)
    return exons


def introns(exons: Dict[str, List[Interval]]) -> Dict[str, List[Interval]]:
    introns = defaultdict(list)
    for tid, regions in exons.items():
        assert len(regions) >= 1
        if len(regions) == 1:
            continue

        regions = sorted(regions, key=lambda x: x.start)
        for i in range(len(regions) - 1):
            assert regions[i].chrom == regions[i + 1].chrom and \
                   regions[i].end <= regions[i + 1].start and \
                   regions[i].strand == regions[i + 1].strand
            if regions[i].end != regions[i + 1].start:
                introns[tid].append(Interval(regions[i].chrom, regions[i].end, regions[i + 1].start, name=tid,
                                             strand=regions[i].strand))
    return introns


def edited_repeats(gnames: Set[str], location=("utr3", "exons", "utr5")) -> Set[Tuple[str, str]]:
    df = pd.read_csv(utils.paths.RESULTS.joinpath("editing-distribution.csv"))
    df = df[
        df['Location'].isin(location) &
        (df['Gene'].apply(lambda x: any(g in gnames for g in x.split(';'))))
        ]
    EDITED = set()
    for coord, cls in zip(df['Coordinates'], df['RepeatCls']):
        chrom, buf = coord.split(":")
        start, end = buf.split('-')
        EDITED.add(utils.miscellaneous.hashinter(Interval(chrom, int(start), int(end)), cls))

    # There are several fused repeats in the EditingIndexer annotation
    # untangle them manually
    fused = {
        'chr17:46774807-46775104': ("chr17:46774807-46774958", "chr17:46774958-46775104"),
        'chr2:121492741-121492924': ("chr2:121492741-121492814", "chr2:121492783-121492924"),
        'chr4:155621389-155621706': (
        "chr4:155621389-155621403", "chr4:155621403-155621550", "chr4:155621550-155621706"),
        'chr5:142870026-142870298': ("chr5:142870026-142870150", "chr5:142870150-142870298"),
        "chr17:35858473-35858692": ("chr17:35858473-35858565", "chr17:35858565-35858692")
    }
    for f, reps in fused.items():
        f = (f, "SINE")
        if f in EDITED:
            EDITED.remove(f)
            EDITED.update({(r, "SINE") for r in reps})

    return EDITED


def edited_genes(location=("utr3", "exons", "utr5")) -> Set[str]:
    df = pd.read_csv(utils.paths.RESULTS.joinpath("editing-distribution.csv"))
    df = df[df['Location'].isin(location)]
    genes = set(chain(*[x.split(';') for x in df['Gene'] if x != "na"]))
    return genes


def expressed_genes(threshold: float = 1e-6, gtype: Optional[str] = None) -> Set[str]:
    expressed = pd.read_csv(paths.DE_GENES.joinpath('expression.normalized.csv'), index_col=0).T.reset_index()
    expressed['IFNb'] = expressed['index'].apply(lambda x: x.split('_')[2])
    expressed = expressed.drop(columns=['index'])

    # mask = expressed[expressed['IFNb'] == f"IFNb-{ifnb}"].median(axis=0) > threshold
    mask = (expressed.groupby(['IFNb']).max() > threshold).any()

    if gtype:
        genes = []
        for gid in mask.index[mask]:
            try:
                if utils.ensembl.gene_id_to_type(gid) == gtype:
                    genes.append(gid)
            except ValueError:
                pass
    else:
        genes = set(mask.index[mask])
    return set(genes)


def repeats_in(regions: Dict[str, List[Interval]]) -> Dict[str, Set[Interval]]:
    index = _repeats_index()
    results = defaultdict(set)
    for tid, intervals in regions.items():
        # map each interval to the repeats index
        for feature in intervals:
            coord = GenomicInterval(feature.chrom, feature.start, feature.end)
            # gather mapped repeats
            for iv, intersected in index[coord].steps():
                if iv.length <= 0:
                    continue
                for repeat in intersected:
                    repname, family, cls = utils.repmasker.classify(repeat.name)
                    # Map low complexity / simple repeats to the region strand
                    assert repeat.strand != "." and feature.strand != "."
                    if repeat.strand != feature.strand:
                        if repeat.name in {"polypyrimidine", "polypurine"}:
                            repeat.name = {"polypyrimidine": "polypurine",
                                           "polypurine": "polypyrimidine"}[repeat.name]
                            repeat.strand = feature.strand
                        elif cls == "Low_complexity":
                            if repname.endswith("-rich"):
                                suffix = "-rich"
                            else:
                                assert repname.endswith("_rich"), repname
                                suffix = "_rich"
                            repname = repname.replace(suffix, "")
                            repname = Seq(repname).reverse_complement()
                            repname = f"{repname}{suffix}"
                            repeat.name = repname
                            repeat.strand = feature.strand
                        elif cls == "Simple_repeat":
                            assert repname[0] == "(" and repname[-2:] == ")n", repname
                            repname = Seq(repname[1:-2]).reverse_complement()
                            repname = f"({repname})n"
                            repeat.name = repname
                            repeat.strand = feature.strand
                    # Save the result
                    results[tid].add(repeat)
    return results


def match_inverted(
        repeats: List[Interval], prioritized: Set[Tuple[str, str]], disthr: int = 3_000
) -> Tuple[List[Tuple[Interval, Interval]], List[Interval]]:
    if len(set(x.strand for x in repeats)) == 1:
        return [], repeats
    result, left = [], []

    minus = [x for x in repeats if x.strand == "-"]
    plus = [x for x in repeats if x.strand == "+"]

    db, query = (minus, plus) if len(minus) > len(plus) else (plus, minus)
    # Presort to match edited repeats first
    query = sorted(query, key=lambda x: hashinter(x, "SINE") in prioritized, reverse=True)
    for r in query:
        # Match edited first, then by distance
        suitable = sorted(
            db,
            key=lambda x: (0 if hashinter(x, "SINE") in prioritized else 1, distance(x, r))
        )[0]
        if distance(suitable, r) < disthr:
            db.remove(suitable)
            result.append((r, suitable))
        else:
            left.append(r)
    left.extend(db)
    return result, left
