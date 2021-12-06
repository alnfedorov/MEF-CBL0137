from pathlib import Path

RESOURCES = Path(__file__).parent.parent.parent.joinpath("resources").absolute()
RESULTS = Path(__file__).parent.parent.parent.joinpath("results").absolute()

SALMON = RESOURCES.joinpath("Salmon")
SIGNAL = RESOURCES.joinpath("signal")

DE_GENES = RESULTS.joinpath("DE-genes")
DE_CLUSTERS = RESULTS.joinpath("DE-clusters")

REPMASKER = RESOURCES.joinpath("mm10-repeat-masker.bed.gz")

EDITING_INDEX_ROOT = RESOURCES.joinpath("editing-index")
EDITING_INDEX_INTERVALS = EDITING_INDEX_ROOT.joinpath("intervals")
EDITING_INDEX_VALUES = EDITING_INDEX_ROOT.joinpath("values")

CURATED_EDITED_REPEATS = RESOURCES.joinpath("edited-repeats.txt")
CURATED_MAYBE_EDITED_REPEATS = RESOURCES.joinpath("maybe-edited-repeats.txt")
CURATED_NOT_EDITED_REPEATS = RESOURCES.joinpath("not-edited-repeats.txt")

GENCODE_5UTR = RESOURCES.joinpath("5`utr-gencode-vm25.bed.gz")
GENCODE_EXONS = RESOURCES.joinpath("exons-gencode-vm25.bed.gz")
GENCODE_INTRONS = RESOURCES.joinpath("introns-gencode-vm25.bed.gz")
GENCODE_3UTR = RESOURCES.joinpath("3`utr-gencode-vm25.bed.gz")
GENCODE_GENES = RESOURCES.joinpath("genes-gencode-vm25.bed.gz")

FRAGMENTS_TO_FEATURES = RESULTS.joinpath("fragments-to-exons-introns-intergenic")
FRAGMENTS_TO_GENES = RESULTS.joinpath("fragments-to-genes")

BAM = RESOURCES.joinpath("BAM")