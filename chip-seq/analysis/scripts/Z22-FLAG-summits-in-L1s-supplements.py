import pathlib
from typing import Tuple

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
        if (enrichment[:utr5end].max() / enrichment.max()) < 0.75:
            utr5_enriched = False
            break
    mask.append(utr5_enriched)
mask = np.asarray(mask)
print(f"5`UTR enriched: {sum(mask)}/{mask.size}")
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
    conssummit = utils.miscellaneous.repeat_position_to_consensus_position(matching, seqsummit)
    return conssummit, mm10summit


conssummits, mm10summits = [], []
for _, row in df.iterrows():
    conssummit, mm10summit = infer_summit(
        row['repeat'], [z22enrich, flagenrich], row['matching-to-consensus'], row['chrom'],
        row['repeat-start(mm10)'], row['repeat-end(mm10)'], row['strand']
    )
    conssummits.append(conssummit)
    mm10summits.append(mm10summit)

df['Z22-FLAG-enrichment-5`UTR-summit(consensus)'] = conssummits
df['Z22-FLAG-enrichment-5`UTR-summit(mm10)'] = mm10summits

# 5 - get a roi for each summit
starts, ends = [], []
offset = 150
for _, row in df.iterrows():
    summit = row['Z22-FLAG-enrichment-5`UTR-summit(mm10)']
    starts.append(
        max(row['repeat-start(mm10)'], summit - offset)
    )
    ends.append(
        min(row['repeat-end(mm10)'], summit + offset)
    )

df['5`UTR-summit-window-start(mm10)'] = starts
df['5`UTR-summit-window-end(mm10)'] = ends


# 5 - get sequence for each window
windows = []
for _, row in df.iterrows():
    windows.append(
        Interval(row['chrom'], row['5`UTR-summit-window-start(mm10)'], row['5`UTR-summit-window-end(mm10)'], strand=row['strand'])
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

df['5`UTR-summit-window-seq'] = window_seq

# 5 - run zhunt for each window
bestzseq, bestzscore, bestzpp, bestzseqstart, bestzseqend = [], [], [], [], []
for sequence in df['5`UTR-summit-window-seq']:
    result = utils.zhunt.run(sequence)
    ind = result['Z-Score'].argmax()

    bestzseq.append(result['Sequence'][ind].upper())
    bestzscore.append(result['Z-Score'][ind])
    bestzpp.append(result['Conformation'][ind])

    start, end = result['Start'][ind] - 1, result['End'][ind] - 1
    bestzseqstart.append(start)
    bestzseqend.append(end)

df['best-zseq-in-window'] = bestzseq
df['best-zseq-in-window-zscore'] = bestzscore
df['best-zseq-in-window-conformation'] = bestzpp
df['best-zseq-in-window-start(window)'] = bestzseqstart
df['best-zseq-in-window-end(window)'] = bestzseqend

# 6 - map zhunt window to the consensus
conscenter, mm10start, mm10end = [], [], []
for _, row in df.iterrows():
    window_start, window_end = row['5`UTR-summit-window-start(mm10)'], row['5`UTR-summit-window-end(mm10)']
    repeat_start, repeat_end = row['repeat-start(mm10)'], row['repeat-end(mm10)']
    if row['strand'] == "+":
        offset = window_start - repeat_start
    else:
        assert row['strand'] == "-"
        offset = repeat_end - window_end

    assert offset >= 0
    # repeat coordinates
    zstart, zend = row['best-zseq-in-window-start(window)'], row['best-zseq-in-window-end(window)']
    zstart, zend = zstart + offset, zend + offset
    # consensus coordinates
    consstart = utils.miscellaneous.repeat_position_to_consensus_position(row['matching-to-consensus'], zstart)
    consend = utils.miscellaneous.repeat_position_to_consensus_position(row['matching-to-consensus'], zend)
    conscenter.append((consstart + consend) / 2)

    # mm10 coordinates
    if row['strand'] == "+":
        mm10start.append(repeat_start + zstart)
        mm10end.append(repeat_start + zend)
    else:
        mm10start.append(repeat_end - zend)
        mm10end.append(repeat_end - zstart)

df['best-zseq-in-window-center(consensus)'] = conscenter
df['best-zseq-in-window-start(mm10)'] = mm10start
df['best-zseq-in-window-end(mm10)'] = mm10end


# 7 - get best z-score for the whole 5`UTR
bestzseq, bestzscore, bestzpp, conscenter, mm10start, mm10end = [], [], [], [], [], []
for _, row in df.iterrows():
    result = utils.zhunt.run(row['repeat-sequence'])
    ind = result['Z-Score'].argmax()

    bestzseq.append(result['Sequence'][ind].upper())
    bestzscore.append(result['Z-Score'][ind])
    bestzpp.append(result['Conformation'][ind])

    zstart, zend = result['Start'][ind] - 1, result['End'][ind] - 1
    consstart = utils.miscellaneous.repeat_position_to_consensus_position(row['matching-to-consensus'], zstart)
    consend = utils.miscellaneous.repeat_position_to_consensus_position(row['matching-to-consensus'], zend)
    conscenter.append((consstart + consend) / 2)

    # mm10 coordinates
    if row['strand'] == "+":
        mm10start.append(row['repeat-start(mm10)'] + zstart)
        mm10end.append(row['repeat-start(mm10)'] + zend)
    else:
        mm10start.append(row['repeat-end(mm10)'] - zend)
        mm10end.append(row['repeat-end(mm10)'] - zstart)

df['best-zseq-in-repeat'] = bestzseq
df['best-zseq-in-repeat-zscore'] = bestzscore
df['best-zseq-in-repeat-conformation'] = bestzpp
df['best-zseq-in-repeat-start(mm10)'] = mm10start
df['best-zseq-in-repeat-end(mm10)'] = mm10end
df['best-zseq-in-repeat-center(consensus)'] = conscenter

df = df[[
    'repeat', 'chrom', 'repeat-start(mm10)', 'repeat-end(mm10)', 'strand',
    'Z22-FLAG-enrichment-5`UTR-summit(consensus)', '5`UTR-summit-window-start(mm10)', '5`UTR-summit-window-end(mm10)',
    # summit window z-sequence
    '5`UTR-summit-window-seq', 'best-zseq-in-window', 'best-zseq-in-window-zscore', 'best-zseq-in-window-conformation',
    'best-zseq-in-window-start(mm10)', 'best-zseq-in-window-end(mm10)', 'best-zseq-in-window-center(consensus)',
    # repeat window z-sequence
    'best-zseq-in-repeat', 'best-zseq-in-repeat-zscore', 'best-zseq-in-repeat-conformation',
    'best-zseq-in-repeat-start(mm10)', 'best-zseq-in-repeat-end(mm10)', 'best-zseq-in-repeat-center(consensus)',
]]

saveto = pathlib.Path(__file__).name.replace(".py", ".tsv")
saveto = utils.paths.RESULTS.joinpath(saveto)
df.to_csv(saveto, sep="\t", index=False)
