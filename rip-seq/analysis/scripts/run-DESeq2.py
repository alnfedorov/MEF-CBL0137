import os
from pathlib import Path
from subprocess import check_call

import pandas as pd

import utils

# Prefilter Salmon results to remove rRNA/tRNA genes
# It might be not needed, but it won't hurt to ensure that they don't influence the process at all
toignore = ["ENSMUSG00000064337", "ENSMUSG00000064339", "ENSMUSG00000098178", "ENSMUSG00000035202",
            "ENSMUSG00000106106", "ENSMUSG00000099250", "ENSMUSG00000065037", "ENSMUSG00000097971",
            "ENSMUSG00000099021"]
for file in utils.paths.SALMON.iterdir():
    if not file.name.endswith(".quant.genes.sf"):
        continue
    df = pd.read_csv(file, sep="\t").set_index("Name")
    mask = [x not in toignore for x in df.index]
    df = df[mask]
    df.to_csv(file, sep="\t")

utils.paths.DE_GENES.mkdir(parents=True, exist_ok=True)
script = Path(__file__).absolute().as_posix().replace(".py", ".R")
check_call(["Rscript", script, utils.paths.SALMON, utils.paths.DE_GENES])

for file in os.listdir(utils.paths.DE_GENES):
    file = utils.paths.DE_GENES.joinpath(file)
    df = pd.read_csv(file, index_col=0)

    df.index.rename("Ensembl gene ID", inplace=True)
    names = []
    for geneid in df.index:
        try:
            names.append(utils.ensembl.id_to_name(geneid))
        except KeyError:
            names.append(None)
    df['Gene name'] = names
    df.to_csv(file)
