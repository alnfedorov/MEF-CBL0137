import numpy as np
import parasail
from pybedtools import BedTool

from . import paths


def is_valid_chromosome(chrom: str):
    return not ("random" in chrom or "Un" in chrom)


def match_to_consensus(consensus: str, repeat: str, gapopen: int = 5, gapext: int = 2):
    """
    Calculate matching: consensus position -> repeat position (-1 if matching does not exist = deletion in the repeat)
    """
    matrix = parasail.Matrix(paths.SUBSTITUTION_MATRIX.as_posix(), case_sensitive=True)

    template = consensus.upper()
    profile = parasail.profile_create_16(template, matrix)

    # Align to the consensus
    query = repeat.upper()
    cigar = parasail.nw_trace_scan_profile_16(profile, query, gapopen, gapext).cigar

    qprev, tprev, step = 0, 0, []
    matching = np.full(len(template), -1, dtype=np.int32)
    for cigar_char in cigar.decode.decode():
        if cigar_char.isdigit():  # still number
            step.append(cigar_char)
        else:  # cigar op
            step = int(''.join(step))
            if cigar_char == "=":  # match
                matching[tprev: tprev + step] = np.arange(qprev, qprev + step)
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
    return matching


def repeat_position_to_consensus_position(matching: np.ndarray, repind: int) -> int:
    """Find closest matched consensus position for the given repeat coordinate"""
    ind = np.abs(matching - repind)
    ind[matching == -1] = ind.max()
    ind = ind.argmin()
    assert matching[ind] != -1
    return ind


def annotate_rectangles_with_values(rectangles, ax):
    # https://matplotlib.org/stable/gallery/lines_bars_and_markers/barchart.html
    for rect in rectangles:
        height = rect.get_height()
        ax.annotate(f"{height}%",
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 1),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
