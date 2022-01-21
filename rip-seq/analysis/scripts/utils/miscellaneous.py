from typing import List, Tuple

from pybedtools import Interval


def annotate_rectangles_with_values(rectangles, ax):
    # https://matplotlib.org/stable/gallery/lines_bars_and_markers/barchart.html
    for rect in rectangles:
        height = rect.get_height()
        ax.annotate(f"{height}%",
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 1),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')


def merge(input: List[Interval]):
    if len(input) == 0:
        return []
    input.sort(key=lambda x: x.start)

    # array to hold the merged intervals
    merged = []
    chrom, start, end = input[0].chrom, input[0].start, input[0].end
    for current in input:
        if current.start > end:
            merged.append(Interval(chrom, start, end))
            end = current.end
            start = current.start
        elif current.end > end:
            end = current.end
    merged.append(Interval(chrom, start, end))
    return merged


def distance(x: Interval, y: Interval) -> float:
    return abs((x.start + x.end) / 2 - (y.start + y.end) / 2)


# Custom __hash__ for Interval (pybedtools)
def hashinter(x: Interval, cls: str) -> Tuple[str, str]:
    return (f"{x.chrom}:{x.start}-{x.end}", cls)
