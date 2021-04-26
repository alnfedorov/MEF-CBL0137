from pathlib import Path

RESOURCES_FOLDER = Path(__file__).parent.parent.parent.joinpath("resources").absolute()

MM10_GENOME = RESOURCES_FOLDER.joinpath("genome.fa")

REPEATMASKER = RESOURCES_FOLDER.joinpath("mm10-repeat-masker.bed.gz")
BLACKLIST = RESOURCES_FOLDER.joinpath("mm10-blacklist.v2.bed.gz")

L1BASE_FULL_LEN_NON_INTACT = RESOURCES_FOLDER.joinpath("L1Base/full-length-non-intact.bed")
L1BASE_ORF2_INTACT = RESOURCES_FOLDER.joinpath("L1Base/ORF2-intact.bed")
L1BASE_ORF1_ORF2_INTACT = RESOURCES_FOLDER.joinpath("L1Base/ORF1&ORF2-intact.bed")

Z22_CURAX_14H = RESOURCES_FOLDER.joinpath("peaks/curax-14h-Z22-vs-IgG-macs2.narrowPeak")
Z22_CURAX_14H_FE = RESOURCES_FOLDER.joinpath("signal/curax-14h-Z22-vs-IgG-macs2-fe.bw")

FLAG_CURAX_0H = RESOURCES_FOLDER.joinpath("peaks/curax-0h-FLAG-vs-IgG-macs2.narrowPeak")
FLAG_CURAX_14H = RESOURCES_FOLDER.joinpath("peaks/curax-14h-FLAG-vs-IgG-macs2.narrowPeak")
FLAG_CURAX_14H_FE = RESOURCES_FOLDER.joinpath("signal/curax-14h-FLAG-vs-IgG-macs2-fe.bw")

SUBSTITUTION_MATRIX = RESOURCES_FOLDER.joinpath("substitution-matrix.txt")

RESULTS = Path(__file__).parent.parent.parent.joinpath("results").absolute()
