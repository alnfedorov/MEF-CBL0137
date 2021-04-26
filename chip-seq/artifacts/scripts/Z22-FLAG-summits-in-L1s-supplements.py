import numpy as np
from tqdm import tqdm
from pybedtools import BedTool, Interval
from Bio import SeqIO
import pyBigWig
import pandas as pd

import utils


# 1 - fetch all relevant repeats
blacklisted = BedTool(utils.paths.BLACKLIST)

z22 = BedTool(utils.paths.Z22_CURAX_14H).subtract(blacklisted, A=True)
flag = BedTool(utils.paths.FLAG_CURAX_14H).subtract(blacklisted, A=True)

z22inter = z22.intersect(flag, u=True, wa=True)
flaginter = flag.intersect(z22, u=True, wa=True)
roi = z22inter.cat(flaginter, postmerge=True)

intactl1s = BedTool(utils.paths.L1BASE_ORF1_ORF2_INTACT).sort()
repmasker = BedTool(utils.paths.REPEATMASKER).sort()\
    .intersect(intactl1s, u=True, wa=True)\
    .intersect(roi, u=True, wa=True)

REPEATS = {"L1Md_A": [], "L1Md_T": []}
for r in tqdm(repmasker):
    # keep only highly conserved sequences
    if r.name in REPEATS and 6000 <= r.length <= 8000:
        REPEATS[r.name].append(r)
print({k: len(v) for k, v in REPEATS.items()})
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
                   "repeat-start": starts, "repeat-end": ends, "strand": strands})

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
        enrichment = bigwig.values(row['chrom'], row['repeat-start'], row['repeat-end'], numpy=True)
        if row['strand'] == "-":
            enrichment = enrichment[::-1]

        # map to consensus coordinates
        enrichment = np.asarray([enrichment[x] if x != -1 else 0 for x in row['matching-to-consensus']])

        utr5end = utils.l1consensus.structure(row['repeat'])["5`UTR"][1]
        if (enrichment[:utr5end].max() / enrichment.max()) < 0.75:
            utr5_enriched = False
            break
    mask.append(utr5_enriched)
mask = np.asarray(mask)
print(f"5`UTR enriched: {sum(mask)}/{mask.size}")
df = df[mask]


# 5 - get joint Z22 and FLAG enrichment summit for each repeat
def infer_summit(repeat: str, bigwigs: list, matching: np.ndarray,
                 chrom: str, start: int, end: int, strand: str) -> int:
    utr5end = utils.l1consensus.structure(repeat)["5`UTR"][1]
    utr5 = matching[:utr5end].max()    # end of the 5`UTR region in the repeat coordinates
    assert utr5 >= 0
    if strand == "+":
        values = [b.values(chrom, start, start + utr5, numpy=True) for b in bigwigs]
        # L1-normalize FLAG and Z22 enrichment
        values = [x / x.sum() for x in values]
        # Pick the region where joint enrichment is the highest
        values = np.sum(values, axis=0)
        peak = values.argmax()
    else:
        assert strand == "-"
        values = [b.values(chrom, end - utr5, end, numpy=True) for b in bigwigs]
        # Same idea as above
        values = [x / x.sum() for x in values]
        values = np.sum(values, axis=0)
        peak = values.argmax()
        peak = utr5 - peak

    peak = utils.miscellaneous.repeat_position_to_consensus_position(matching, peak)
    return peak


summits = []
for _, row in df.iterrows():
    summits.append(
        infer_summit(row['repeat'], [z22enrich, flagenrich], row['matching-to-consensus'], row['chrom'],
                     row['repeat-start'], row['repeat-end'], row['strand'])
    )
df['Z22-FLAG-enrichment-summit-in-consensus'] = summits

# 5 - get a roi for each summit
starts, ends = [], []
offset = 150
for _, row in df.iterrows():
    summit = row['Z22-FLAG-enrichment-summit-in-consensus']
    window = max(0, summit - offset), summit + offset
    window = row['matching-to-consensus'][window[0]: window[1]]
    window = window[window >= 0]
    start, end = window.min(), window.max()

    if row['strand'] == "+":
        starts.append(row['repeat-start'] + start)
        ends.append(row['repeat-start'] + end)
    else:
        starts.append(row['repeat-end'] - end)
        ends.append(row['repeat-end'] - start)

    assert starts[-1] < ends[-1] and \
           row['repeat-start'] <= starts[-1] <= row['repeat-end'] and \
           row['repeat-start'] <= ends[-1] <= row['repeat-end']

df['summit-window-start'] = starts
df['summit-window-end'] = ends

# 5 - get sequence for each window
windows = []
for _, row in df.iterrows():
    windows.append(
        Interval(row['chrom'], row['summit-window-start'], row['summit-window-end'], strand=row['strand'])
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

    window_seq.append(str(seq.seq))

df['summit-window-sequence'] = window_seq

# 5 - run zhunt for each window
bestzseq, bestzscore, bestzpp, bestzseqstart, bestzseqend = [], [], [], [], []
for sequence in df['summit-window-sequence']:
    result = utils.zhunt.run(sequence)
    ind = result['Z-Score'].argmax()

    bestzseq.append(result['Sequence'][ind].upper())
    bestzscore.append(result['Z-Score'][ind])
    bestzpp.append(result['Conformation'][ind])

    start, end = result['Start'][ind] - 1, result['End'][ind] - 1
    bestzseqstart.append(start)
    bestzseqend.append(end)

df['best-zseq-in-window'] = bestzseq
df['best-zscore-in-window'] = bestzscore
df['best-zconformation'] = bestzpp
df['best-zseq-in-window-start'] = bestzseqstart
df['best-zseq-in-window-end'] = bestzseqend

# 6 - map zhunt window to the consensus
consstart, consend = [], []
for _, row in df.iterrows():
    window_start, window_end = row['summit-window-start'], row['summit-window-end']
    repeat_start, repeat_end = row['repeat-start'], row['repeat-end']
    if row['strand'] == "+":
        offset = window_start - repeat_start
    else:
        assert row['strand'] == "-"
        offset = repeat_end - window_end

    assert offset >= 0
    # repeat coordinates
    zstart, zend = row['best-zseq-in-window-start'], row['best-zseq-in-window-end']
    zstart, zend = zstart + offset, zend + offset

    consstart.append(
        utils.miscellaneous.repeat_position_to_consensus_position(row['matching-to-consensus'], zstart)
    )
    consend.append(
        utils.miscellaneous.repeat_position_to_consensus_position(row['matching-to-consensus'], zend)
    )

df['best-zseq-in-consensus-start'] = consstart
df['best-zseq-in-consensus-end'] = consend
