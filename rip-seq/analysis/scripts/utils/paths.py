from pathlib import Path

RESOURCES = Path(__file__).parent.parent.parent.joinpath("resources").absolute()
RESULTS = Path(__file__).parent.parent.parent.joinpath("results").absolute()

INTERFEROME_QUERIES = RESOURCES.joinpath("Interferome")
ISG_INTERFEROME = RESULTS.joinpath("interferome-all.csv")

ISG_MEF_MICROARRAY = RESOURCES.joinpath("IFNb.treatment.2fold (ISGs list).xlsx")
ISG_MEF_NGS = RESOURCES.joinpath("IFNb-MEF-NGS.csv")

SALMON = RESOURCES.joinpath("Salmon")
SIGNAL = RESOURCES.joinpath("signal")

DE_GENES = RESULTS.joinpath("DE-genes")
ADAR1_WT_IFNb_72h_Z22_vs_IgG = DE_GENES.joinpath("ADAR1-WT-72h-Z22-vs-IgG.csv")

REPMASKER = RESOURCES.joinpath("mm10-repeat-masker.bed.gz")
COVERAGE_BLACKLIST = RESULTS.joinpath("coverage-blacklist.bed.gz")

EDITING_INDEX_ROOT = RESULTS.joinpath("editing-index")
EDITING_INDEX_INTERVALS = EDITING_INDEX_ROOT.joinpath("intervals")
EDITING_INDEX_VALUES = EDITING_INDEX_ROOT.joinpath("values")

GENCODE_5UTR = RESOURCES.joinpath("5`utr-gencode-vm25.bed.gz")
GENCODE_EXONS = RESOURCES.joinpath("exons-gencode-vm25.bed.gz")
GENCODE_INTRONS = RESOURCES.joinpath("introns-gencode-vm25.bed.gz")
GENCODE_3UTR = RESOURCES.joinpath("3`utr-gencode-vm25.bed.gz")
FRAGMENTS_TO_FEATURES = RESULTS.joinpath("fragments-to-exons-introns-intergenic")
FRAGMENTS_TO_GENES = RESULTS.joinpath("fragments-to-genes")
BAM = RESOURCES.joinpath("BAM")
