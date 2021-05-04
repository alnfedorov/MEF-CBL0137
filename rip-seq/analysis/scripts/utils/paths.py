from pathlib import Path

RESOURCES = Path(__file__).parent.parent.parent.joinpath("resources").absolute()
RESULTS = Path(__file__).parent.parent.parent.joinpath("results").absolute()

INTERFEROME_ALL_QUERIES = RESOURCES.joinpath("Interferome/all-genes")
INTERFEROME_ALL = RESULTS.joinpath("interferome-merged.csv")

INTERFEROME_MEF_QUERY = RESOURCES.joinpath("Interferome/mef/results.txt")
INTERFEROME_MEF = RESULTS.joinpath("interferome-mef.csv")

ISG_IN_HOME_MICROARRAYS = RESOURCES.joinpath("IFNb.treatment.2fold (ISGs list).xlsx")

SALMON = RESOURCES.joinpath("Salmon")
DE_GENES = RESULTS.joinpath("DE-genes")
ADAR1_WT_IFNb_72h_Z22_vs_IgG = DE_GENES.joinpath("ADAR1-WT-72h-Z22-vs-IgG.csv")
