import pathlib
from typing import Tuple
from collections import Counter

import numpy as np
import pandas as pd
import pyBigWig
from Bio import SeqIO
from pybedtools import BedTool, Interval
from tqdm import tqdm

import utils

# 1 - fetch all relevant repeats
blacklisted = BedTool(utils.paths.BLACKLIST)

z22 = BedTool(utils.paths.Z22_CURAX_14H).subtract(blacklisted, A=True)
flag = BedTool(utils.paths.FLAG_CURAX_14H).subtract(blacklisted, A=True)

z22inter = z22.intersect(flag, u=True, wa=True)
flaginter = flag.intersect(z22, u=True, wa=True)
roi = z22inter.cat(flaginter, postmerge=True)

intactl1s = BedTool(utils.paths.L1BASE_ORF1_ORF2_INTACT).sort()
repmasker = BedTool(utils.paths.REPEATMASKER).sort() \
    .intersect(intactl1s, u=True, wa=True) \
    .intersect(roi, u=True, wa=True)

REPEATS = {"L1Md_A": [], "L1Md_T": []}
for r in tqdm(repmasker):
    # keep only highly conserved sequences
    if r.name in REPEATS:
        REPEATS[r.name].append(r)
print("Total repeats: ", {k: len(v) for k, v in REPEATS.items()})
REPEATS = {k: BedTool(v).sort() for k, v in REPEATS.items()}

# 2 - fetch sequences for each repeat and create a dataframe
reptypes, sequences, chroms, starts, ends, strands = [], [], [], [], [], []
for repeat, bed in REPEATS.items():
    fasta = bed.sequence(
        fi=utils.paths.MM10_GENOME.as_posix(), s=True
    ).seqfn
    for seq in SeqIO.parse(fasta, format="fasta"):
        buffer, strand = seq.id[:-3], seq.id[-3:]
        chrom, buffer = buffer.split(":")
        start, end = buffer.split("-")
        assert strand in ("(-)", "(+)")
        strand = strand[1]
        reptypes.append(repeat)
        sequences.append(str(seq.seq))
        chroms.append(chrom)
        starts.append(int(start))
        ends.append(int(end))
        strands.append(strand)

df = pd.DataFrame({"repeat": reptypes, "repeat-sequence": sequences, "chrom": chroms,
                   "repeat-start(mm10)": starts, "repeat-end(mm10)": ends, "strand": strands})

# 3 - create a mapping for each repeat from the consensus coordinates to repeat coordinates
matchings = []
for sequence, repeat in zip(df['repeat-sequence'], df['repeat']):
    consensus = utils.l1consensus.sequence(repeat)
    matchings.append(
        utils.miscellaneous.match_to_consensus(consensus, sequence)
    )
df['matching-to-consensus'] = matchings

# 4 - remove repeats with enrichment summit outside the 5`UTR
z22enrich = pyBigWig.open(utils.paths.Z22_CURAX_14H_FE.as_posix())
flagenrich = pyBigWig.open(utils.paths.FLAG_CURAX_14H_FE.as_posix())

mask = []
for _, row in df.iterrows():
    utr5_enriched = True
    for bigwig in [z22enrich, flagenrich]:
        enrichment = bigwig.values(row['chrom'], row['repeat-start(mm10)'], row['repeat-end(mm10)'], numpy=True)
        if row['strand'] == "-":
            enrichment = enrichment[::-1]

        # map to consensus coordinates
        enrichment = np.asarray([enrichment[x] if x != -1 else 0 for x in row['matching-to-consensus']])

        utr5end = utils.l1consensus.structure(row['repeat'])["5`UTR"][1]
        if enrichment[:utr5end].max() != enrichment.max():
            utr5_enriched = False
            break
    mask.append(utr5_enriched)
mask = np.asarray(mask)
print(f"5`UTR enriched: {mask.sum()}/{mask.size}")
df = df[mask]


# 5 - get joint Z22 and FLAG enrichment summit for each repeat
def infer_summit(repeat: str, bigwigs: list, matching: np.ndarray,
                 chrom: str, start: int, end: int, strand: str) -> Tuple[int, int]:
    utr5end = utils.l1consensus.structure(repeat)["5`UTR"][1]
    utr5 = matching[:utr5end].max()  # end of the 5`UTR region in the repeat coordinates
    assert utr5 >= 0
    if strand == "+":
        values = [b.values(chrom, start, start + utr5, numpy=True) for b in bigwigs]
        # L1-normalize FLAG and Z22 enrichment
        values = [x / x.sum() for x in values]
        # Pick the region where joint enrichment is the highest
        values = np.sum(values, axis=0)
        seqsummit = values.argmax()
        # to mm10
        mm10summit = start + seqsummit
    else:
        assert strand == "-"
        values = [b.values(chrom, end - utr5, end, numpy=True) for b in bigwigs]
        # Same idea as above
        values = [x / x.sum() for x in values]
        values = np.sum(values, axis=0)[::-1]
        seqsummit = values.argmax()
        mm10summit = end - seqsummit
    return mm10summit


conssummits, mm10summits = [], []
for _, row in df.iterrows():
    mm10summit = infer_summit(
        row['repeat'], [z22enrich, flagenrich], row['matching-to-consensus'], row['chrom'],
        row['repeat-start(mm10)'], row['repeat-end(mm10)'], row['strand']
    )
    mm10summits.append(mm10summit)

df['joint-5`UTR-enrichment-summit(mm10)'] = mm10summits

# 5 - get a roi for each summit
starts, ends = [], []
SUMMIT_OFFSET = 150
for _, row in df.iterrows():
    summit = row['joint-5`UTR-enrichment-summit(mm10)']
    starts.append(
        max(row['repeat-start(mm10)'], summit - SUMMIT_OFFSET)
    )
    ends.append(
        min(row['repeat-end(mm10)'], summit + SUMMIT_OFFSET)
    )

df['joint-summit-150bp-window-start(mm10)'] = starts
df['joint-summit-150bp-window-end(mm10)'] = ends


# 5 - get sequence for each window
windows = []
for _, row in df.iterrows():
    windows.append(
        Interval(row['chrom'], row['joint-summit-150bp-window-start(mm10)'], row['joint-summit-150bp-window-end(mm10)'], strand=row['strand'])
    )

fasta = BedTool(windows).sequence(
    fi=utils.paths.MM10_GENOME.as_posix(), s=True
).seqfn

window_seq = []
for window, seq in zip(windows, SeqIO.parse(fasta, format="fasta")):
    buffer, strand = seq.id[:-3], seq.id[-3:]
    chrom, buffer = buffer.split(":")
    start, end = buffer.split("-")
    assert strand in ("(-)", "(+)")
    strand = strand[1]
    assert window.chrom == chrom and window.start == int(start) and window.end == int(end) and window.strand == strand

    window_seq.append(str(seq.seq).upper())

df['joint-summit-window-seq'] = window_seq


# 5 - run zhunt for each window
def run_zhunt(strand, mm10_repeat_start, mm10_repeat_end, mm10_summit, sequence: str, offset: int = 0) -> dict:
    """
    Run zhunt and return localization and scoring for the best Z-forming sequence.
    In case of multiple sequences, return closest to the enrichment summit.

    offset: distance from the repeat start, >0 for a window inside the repeat sequence
    """
    def tomm10(zstart, zend):
        # repeat coordinates
        zstart, zend = zstart + offset, zend + offset
        # mm10 coordinates
        if strand == "+":
            zstart, zend = mm10_repeat_start + zstart, mm10_repeat_start + zend
        else:
            zstart, zend = mm10_repeat_end - zend, mm10_repeat_end - zstart
        return zstart, zend

    result = utils.zhunt.run(sequence)
    result['Start'] -= 1
    result['Start'] = result['Start'].astype(np.int64)
    result['End'] = result['End'].astype(np.int64)

    # select best Z-windows closest to the enrichment summit
    maxz = result['Z-Score'].max()
    candind = (result['Z-Score'] == maxz).values.nonzero()[0]
    start, end = tomm10(result['Start'][candind].values,
                        result['End'][candind].values)
    center = ((start + end) / 2)

    distance = np.abs(mm10_summit - center)
    assert distance.max() < 15_000, f"{distance}, {mm10_summit}, {center}"
    maxzind = candind[distance.argmin()]
    assert result['Z-Score'][maxzind] == maxz

    start, end = tomm10(result['Start'][maxzind], result['End'][maxzind])
    return {
        "z-score": maxz,
        "z-seq": result['Sequence'][maxzind].upper(),
        "conformation": result['Conformation'][maxzind],
        "mm10-start": start,
        "mm10-end": end
    }


bestzseq, bestzscore, bestzpp, mm10start, mm10end = [], [], [], [], []
for _, row in df.iterrows():
    window_start, window_end = row['joint-summit-150bp-window-start(mm10)'], row['joint-summit-150bp-window-end(mm10)']
    repeat_start, repeat_end = row['repeat-start(mm10)'], row['repeat-end(mm10)']
    if row['strand'] == "+":
        offset = window_start - repeat_start
    else:
        assert row['strand'] == "-"
        offset = repeat_end - window_end

    result = run_zhunt(
        row['strand'], repeat_start, repeat_end, row['joint-5`UTR-enrichment-summit(mm10)'],
        row['joint-summit-window-seq'], offset
    )

    bestzscore.append(result['z-score'])
    bestzseq.append(result['z-seq'])
    bestzpp.append(result['conformation'])
    mm10start.append(result['mm10-start'])
    mm10end.append(result['mm10-end'])

df['best-zseq-in-window'] = bestzseq
df['best-zseq-in-window-zscore'] = bestzscore
df['best-zseq-in-window-conformation'] = bestzpp
df['best-zseq-in-window-start(mm10)'] = mm10start
df['best-zseq-in-window-end(mm10)'] = mm10end


# 7 - get best z-score for the repeat
bestzseq, bestzscore, bestzpp, mm10start, mm10end = [], [], [], [], []
for _, row in df.iterrows():
    result = run_zhunt(
        row['strand'], row['repeat-start(mm10)'], row['repeat-end(mm10)'],
        row['joint-5`UTR-enrichment-summit(mm10)'], row['repeat-sequence']
    )

    bestzscore.append(result['z-score'])
    bestzseq.append(result['z-seq'])
    bestzpp.append(result['conformation'])
    mm10start.append(result['mm10-start'])
    mm10end.append(result['mm10-end'])

df['best-zseq-in-repeat'] = bestzseq
df['best-zseq-in-repeat-zscore'] = bestzscore
df['best-zseq-in-repeat-conformation'] = bestzpp
df['best-zseq-in-repeat-start(mm10)'] = mm10start
df['best-zseq-in-repeat-end(mm10)'] = mm10end

# Best z-sequence within the OFFSET bp from the summit
matched = df['best-zseq-in-window-zscore'] == df['best-zseq-in-repeat-zscore']
print(f"Best Z-score within the {SUMMIT_OFFSET}bp from the enrichment summit: {matched.mean()*100:.2f}%")

# Map best Z-formers to Z22 peaks
zseqs = []
for _, row in df.iterrows():
    zseqs.append(
        Interval(row['chrom'], row['best-zseq-in-repeat-start(mm10)'], row['best-zseq-in-repeat-end(mm10)'])
    )
zseqs = BedTool(zseqs).sort().intersect(roi, wa=True, u=True)
print(f'{len(zseqs)}/{len(df)}({len(zseqs)/len(df)*100:.2f}%) '
      f'best z-formers are located inside the enrichment peak in 5`UTR')

# Count overrepresented sequences
for repeat, group in df.groupby(['repeat']):
    mask = group['best-zseq-in-window-zscore'] == group['best-zseq-in-repeat-zscore']
    group = group[mask]
    print(repeat)
    for mostcommon, count in Counter(group['best-zseq-in-window']).most_common(2):
        print(f"\t{mostcommon} -> {count / len(group) * 100:.2f}%")

df = df[[
    'repeat', 'chrom', 'repeat-start(mm10)', 'repeat-end(mm10)', 'strand',
    'joint-5`UTR-enrichment-summit(mm10)', 'joint-summit-150bp-window-start(mm10)', 'joint-summit-150bp-window-end(mm10)',
    # summit window z-sequence
    'joint-summit-window-seq', 'best-zseq-in-window', 'best-zseq-in-window-zscore', 'best-zseq-in-window-conformation',
    'best-zseq-in-window-start(mm10)', 'best-zseq-in-window-end(mm10)',
    # repeat window z-sequence
    'best-zseq-in-repeat', 'best-zseq-in-repeat-zscore', 'best-zseq-in-repeat-conformation',
    'best-zseq-in-repeat-start(mm10)', 'best-zseq-in-repeat-end(mm10)',
]]

saveto = pathlib.Path(__file__).name.replace(".py", ".tsv")
saveto = utils.paths.RESULTS.joinpath(saveto)
df.to_csv(saveto, sep="\t", index=False)
