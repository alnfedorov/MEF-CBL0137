import pickle
from collections import defaultdict
from itertools import chain
from multiprocessing import Pool, cpu_count
from pathlib import Path

from HTSeq import GenomicInterval, GenomicArrayOfSets, BAM_Reader, pair_SAM_alignments_with_buffer
from pybedtools import BedTool

import utils
from utils import ensembl

utils.paths.FRAGMENTS_TO_FEATURES.mkdir(exist_ok=True, parents=True)

rrna_repeats = BedTool(utils.paths.REPMASKER).filter(lambda x: "LSU" in x.name or "SSU" in x.name).sort()
rrna_repeats = rrna_repeats.slop(b=1000, genome='mm10').sort().merge()

exons = BedTool(utils.paths.GENCODE_EXONS).sort().subtract(rrna_repeats)
introns = BedTool(utils.paths.GENCODE_INTRONS).sort().subtract(exons).subtract(rrna_repeats)

INDEX = GenomicArrayOfSets('auto', stranded=False)
mapchr = lambda x: x.replace("chr", "") if x != "chrM" else "MT"
for name, bed in zip(['intron', 'rRNA'], [introns, rrna_repeats]):
    for interval in bed:
        ginter = GenomicInterval(mapchr(interval.chrom), interval.start, interval.end)
        INDEX[ginter] += name

for interval in exons:
    transcript = interval.name.split("_")[0]
    geneid = ensembl.transcript_to_gene_id(transcript)
    genetype = ensembl.gene_id_to_type(geneid)

    ginter = GenomicInterval(mapchr(interval.chrom), interval.start, interval.end)
    INDEX[ginter] += genetype


def job(bamfile: str):
    print(f"Started => {bamfile}")
    counts = defaultdict(int)
    reader = BAM_Reader(bamfile)
    acceptedop = {"M", "=", "X"}
    for bundle in pair_SAM_alignments_with_buffer(reader, primary_only=True):
        lmate, rmate = bundle
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
    for file in utils.paths.BAM.glob("ADAR1-KO+ZBP1-KO*Z22*.bam"):
        handlers.append(
            pool.apply_async(job, (file.absolute().as_posix(),))
        )
    for h in handlers:
        h.wait()
