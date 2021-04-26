import os
import tempfile
from pybedtools import BedTool
from multiprocessing import cpu_count
from subprocess import check_call
from Bio import SeqIO
import utils
from collections import defaultdict

# gather intact repeats to build the consensus
import utils.l1consensus

repeats = {t: [] for t in ["L1Md_A", "L1Md_T"]}
minlen, maxlen = 6000, 7000

intactl1s = BedTool(utils.paths.L1BASE_ORF1_ORF2_INTACT).sort()
repmasker = BedTool(utils.paths.REPEATMASKER).sort()

for repeat in repmasker.intersect(intactl1s, u=True, wa=True):
    if repeat.name in repeats and maxlen > repeat.length > minlen:
        repeats[repeat.name].append(repeat)

repeats = {k: BedTool(v).sort() for k, v in repeats.items()}

# build the consensus (takes time)
saveto_consensus = {
    "L1Md_A": utils.l1consensus.L1Md_A_CONSENSUS, "L1Md_T": utils.l1consensus.L1Md_T_CONSENSUS
}
for title, l1s in repeats.items():
    print(f"{title} => {len(l1s)}")
    fasta = BedTool(l1s).sort().sequence(
        fi=utils.paths.MM10_GENOME.as_posix(), s=True
    ).seqfn

    fd, saveto_alignment = tempfile.mkstemp()
    os.close(fd)
    check_call([
        "clustalo", "-i", fasta, "--threads", str(cpu_count()), "--verbose", "--force",
        "--seqtype", "dna", "-o", saveto_alignment,
    ])
    os.unlink(fasta)

    alignment = list(SeqIO.parse(saveto_alignment, format="fasta"))
    alignlen = len(alignment[0].seq)
    alignment = {x.id: str(x.seq) for x in alignment}
    os.unlink(saveto_alignment)

    counts = [defaultdict(int) for _ in range(alignlen)]
    for _, sequence in alignment.items():
        for ind, symb in enumerate(sequence):
            counts[ind][symb] += 1
    consensus = [max(position, key=lambda x: position[x] if x != "-" else -1) for position in counts]
    occupancy = [(sum(position.values()) - position['-']) / len(alignment) for position in counts]

    # peak only positions with occupancy > 0.5
    threshold = 0.5
    sequence = [nc for nc, oc in zip(consensus, occupancy) if oc > threshold]
    sequence = "".join(sequence)
    print(f"{title} => {len(sequence)}")

    with open(saveto_consensus[title], 'w') as stream:
        stream.write(sequence)
