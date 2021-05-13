import pickle
from collections import defaultdict
from itertools import chain
from multiprocessing import Pool, cpu_count
from pathlib import Path

import pandas as pd
from HTSeq import GenomicInterval, GenomicArrayOfSets, BAM_Reader, pair_SAM_alignments
from pybedtools import BedTool

import utils

utils.paths.FRAGMENTS_TO_FEATURES.mkdir(exist_ok=True, parents=True)

rrna_repeats = BedTool(utils.paths.REPMASKER).filter(lambda x: "LSU" in x.name or "SSU" in x.name).sort()
rrna_repeats = rrna_repeats.slop(b=1000, genome='mm10').sort().merge()

exons = BedTool(utils.paths.GENCODE_EXONS).sort().subtract(rrna_repeats)
introns = BedTool(utils.paths.GENCODE_INTRONS).sort().subtract(exons).subtract(rrna_repeats)

INDEX = GenomicArrayOfSets('auto', stranded=False)
for name, bed in zip(['exon', 'intron', 'rRNA'], [exons, introns, rrna_repeats]):
    for interval in bed:
        ginter = GenomicInterval(interval.chrom.replace("chr", ""), interval.start, interval.end)
        INDEX[ginter] += name


def job(bamfile: str):
    print(f"Started => {bamfile}")
    counts = defaultdict(int)
    reader = BAM_Reader(bamfile)
    acceptedop = {"M", "=", "X"}
    for bundle in pair_SAM_alignments(reader, bundle=True, primary_only=True):
        if len(bundle) != 1:
            continue
        lmate, rmate = bundle[0]
        if lmate is None or lmate.supplementary or rmate is None or rmate.supplementary:
            continue
        intervals = chain((x.ref_iv for x in lmate.cigar if x.type in acceptedop and x.size > 0),
                          (x.ref_iv for x in rmate.cigar if x.type in acceptedop and x.size > 0))

        features = set()
        for inter in intervals:
            for iv, val in INDEX[inter].steps():
                features |= val

        if len(features) == 0:
            counts['intergenic'] += 1
        elif len(features) > 1:
            counts[tuple(features)] += 1
        else:
            counts[features.pop()] += 1
    filename = Path(bamfile).name.replace(".bam", ".pkl")
    saveto = utils.paths.FRAGMENTS_TO_FEATURES.joinpath(filename)
    with open(saveto, 'wb') as stream:
        pickle.dump(counts, stream)
    print(f"Finished => {saveto}")


with Pool(processes=cpu_count() - 2) as pool:
    handlers = []
    for file in utils.paths.BAM.iterdir():
        handlers.append(
            pool.apply_async(job, (file.absolute().as_posix(),))
        )
    for h in handlers:
        h.wait()

key, exons, introns, intergenic, ambiguous, rRNA = [], [], [], [], [], []
for file in utils.paths.FRAGMENTS_TO_FEATURES.iterdir():
    with open(file.absolute(), 'rb') as stream:
        data = pickle.load(stream)
    key.append(file.name.replace(".namesorted.pkl", ""))
    exons.append(data.pop('exon'))
    introns.append(data.pop('intron'))
    intergenic.append(data.pop('intergenic'))
    rRNA.append(data.pop('rRNA'))
    ambiguous.append(0)
    for k in data:
        assert isinstance(k, tuple) and len(k) >= 2
        if "rRNA" in k:
            continue
        else:
            ambiguous[-1] += data[k]

df = pd.DataFrame({"sample": key, "exons": exons, "introns": introns,
                   "intergenic": intergenic, "ambiguous": ambiguous})

total = df['exons'] + df['introns'] + df['intergenic'] + df['ambiguous']
for k in ['exons', 'introns', 'intergenic', 'ambiguous']:
    df[k] = df[k] / total
df.to_csv(utils.paths.RESULTS.joinpath("fragments-to-exons-introns-intergenic.csv"), index=False)
