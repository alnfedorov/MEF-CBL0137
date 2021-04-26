import matplotlib.pyplot as plt
import numpy as np
import parasail
import pyBigWig
from Bio import SeqIO
from matplotlib import patches as mpatches
from pybedtools import BedTool, Interval
from tqdm import tqdm

import utils

# fetch all relevant repeats
blacklisted = BedTool(utils.paths.BLACKLIST)

z22 = BedTool(utils.paths.Z22_CURAX_14H).subtract(blacklisted, A=True)
flag = BedTool(utils.paths.FLAG_CURAX_14H).subtract(blacklisted, A=True)

z22inter = z22.intersect(flag, u=True, wa=True)
flaginter = flag.intersect(z22, u=True, wa=True)
roi = z22inter.cat(flaginter, postmerge=True)

intactl1s = BedTool(utils.paths.L1BASE_ORF1_ORF2_INTACT).sort()
repmasker = BedTool(utils.paths.REPEATMASKER).sort()
repmasker = repmasker.intersect(intactl1s, u=True, wa=True).intersect(roi, u=True, wa=True)

repeats = {t: [] for t in ["L1Md_A", "L1Md_T"]}
for r in tqdm(repmasker):
    if r.name in repeats:
        repeats[r.name].append(r)
repeats = {k: BedTool(v).sort() for k, v in repeats.items()}

# fetch their sequence
sequences = {k: [] for k in repeats}
for repname, bedtool in repeats.items():
    fasta = bedtool.sort().sequence(
        fi=utils.paths.MM10_GENOME.as_posix(), s=True
    ).seqfn
    for seq in SeqIO.parse(fasta, format="fasta"):
        buffer, strand = seq.id[:-3], seq.id[-3:]
        chrom, buffer = buffer.split(":")
        start, end = buffer.split("-")
        assert strand in ("(-)", "(+)")
        strand = strand[1]
        sequences[repname].append(((chrom, int(start), int(end), strand), str(seq.seq)))


# align to the consensus and aggregate the enrichment
def align(consensus, queries, enrichmentfabrics, gapopen: int = 5, gapext: int = 2):
    matrix = parasail.Matrix(utils.paths.SUBSTITUTION_MATRIX.as_posix(), case_sensitive=False)

    template = consensus.upper()
    profile = parasail.profile_create_16(template, matrix)

    consensusmeta = {k: np.zeros(len(template), dtype=np.float32) for k in enrichmentfabrics}
    matches = np.zeros(len(template), dtype=np.int32)

    for ind, (interval, query) in tqdm(enumerate(queries)):
        # Align to the consensus
        query = str(query).upper()
        cigar = parasail.nw_trace_scan_profile_16(profile, query, gapopen, gapext).cigar

        # Gather enrichment
        interval = Interval(interval[0], interval[1], interval[2], strand=interval[3])
        querymeta = {k: enrichmentfabrics[k](interval) for k in consensusmeta}
        if interval.strand == "-":
            querymeta = {k: v[::-1] for k, v in querymeta.items()}

        # Store matched consensus enrichment
        qprev, tprev, step = 0, 0, []
        for cigar_char in cigar.decode.decode():
            if cigar_char.isdigit():  # still number
                step.append(cigar_char)
            else:  # cigar op
                step = int(''.join(step))
                if cigar_char == "=":  # match
                    matches[tprev: tprev + step] += 1
                    for k in consensusmeta:
                        consensusmeta[k][tprev: tprev + step] += querymeta[k][qprev: qprev + step]
                    tprev += step
                    qprev += step
                elif cigar_char == "X":  # mismatch
                    tprev += step
                    qprev += step
                elif cigar_char == "D":  # deletion
                    qprev += step
                elif cigar_char == "I":  # insertion
                    tprev += step
                else:
                    assert False, f"Unknown operation {cigar_char}"
                step = []
        assert len(step) == 0
    consensusmeta = {k: x / matches.astype(np.float32) for k, x in consensusmeta.items()}
    return consensusmeta, matches


z22enrich = pyBigWig.open(utils.paths.Z22_CURAX_14H_FE.as_posix())
flagenrich = pyBigWig.open(utils.paths.FLAG_CURAX_14H_FE.as_posix())

enrichmentfabrics = {
    'Z22': lambda inter: z22enrich.values(inter.chrom, inter.start, inter.end, numpy=True),
    'FLAG': lambda inter: flagenrich.values(inter.chrom, inter.start, inter.end, numpy=True),
}

for repname, repsequences in sequences.items():
    consensus = utils.l1consensus.sequence(repname)
    consensusmeta, matches = align(consensus, repsequences, enrichmentfabrics)

    consensusmeta = {k: np.asarray(v) for k, v in consensusmeta.items()}

    # smooth values
    window = 8
    smoothed = {}
    for k, metavals in consensusmeta.items():
        smothvals = np.empty_like(metavals)
        for ind, val in enumerate(metavals):
            start, end = max(0, ind - window // 2), min(len(matches) - 1, ind + window // 2)
            weights = matches[start:end] / np.linalg.norm(matches[start:end], ord=1)
            smothvals[ind] = (metavals[start:end] * weights).sum()
        smoothed[k] = smothvals

    # make the plot
    colors = {"Z22": "#FF4474", "FLAG": "#6AD2F4"}
    semitransparent = {"Z22": "#FFDEE6", "FLAG": "#E4F7FD"}
    fontsize = "x-large"
    arrowsy = 3.85
    lw = 1.5
    right_padding = 200

    yticks = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
    xticks = [0, 1000, 2000, 3000, 4000, 5000, 6000]

    fig = plt.figure(figsize=(12, 8))
    ax = fig.gca()

    ax.set_xticks(xticks)
    ax.set_xticklabels([f"{y}" for y in xticks], fontsize='large')

    ax.set_yticks(yticks)
    ax.set_yticklabels([f"{y}" for y in yticks], fontsize='large')

    ax.tick_params(width=lw)

    for k in consensusmeta:
        ax.plot(smoothed[k], color=colors[k])

    ax.axhline(1, color='black', lw=lw)
    start = len(consensus) + 100
    ax.add_patch(
        mpatches.FancyArrowPatch((start, 1), (start, 0), linewidth=lw,
                                 arrowstyle='<->', mutation_scale=20, color='black')
    )
    ax.text(start + 75, 0.5, 'no enrichment', ha='center', va='center', rotation=-90, fontsize='large')

    for struct, (start, end) in utils.l1consensus.structure(repname).items():
        ax.add_patch(
            mpatches.FancyArrowPatch((start, arrowsy), (end, arrowsy), linewidth=lw,
                                     arrowstyle='<->', mutation_scale=20, color='#BBBBBB')
        )
        ax.vlines(end, 0, arrowsy, ls='--', color='#BBBBBB', lw=lw * 1.2)
        ax.vlines(start, 0, arrowsy, ls='--', color='#BBBBBB', lw=lw * 1.2)
        center = (start + end) / 2
        ax.text(center, arrowsy + 0.03, struct, ha='center', va='bottom', fontsize=fontsize, fontstyle='italic')

    ax.set_xlim(0, len(consensus) + right_padding)
    ax.set_ylim((0, yticks[-1] + 0.25))

    ax.set_xlabel("Consensus position", fontsize='xx-large')
    ax.set_ylabel("Mean fold enrichment", fontsize='xx-large')

    legend = [mpatches.Patch(facecolor=colors[k], label=k, edgecolor='black', linewidth=0.75) for k in ['Z22', 'FLAG']]
    ax.legend(bbox_to_anchor=(len(consensus), arrowsy), handles=legend, frameon=False,
              fontsize='xx-large', bbox_transform=ax.transData, loc='upper left')
    ax.set_title(f"Enrichment profile of intact Z22 & FLAG bound $\\bf{{{repname.replace('_', '~')}}}$ repeats",
                 fontsize='xx-large', loc='left')

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    for axis in ['bottom', 'left']:
        ax.spines[axis].set_linewidth(lw)

    fig.tight_layout()
    fig.savefig(utils.paths.RESULTS.joinpath(f"enrichment-profiles-{repname.replace('_', '~')}.svg"),
                bbox_inches="tight", pad_inches=0)
