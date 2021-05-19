import os
from pathlib import Path
from subprocess import check_call

import pandas as pd

import utils

# Lars2 contains LSU rRNA repeat in the 3`UTR.
# It artificially amplifies the signal because Salmon align rRNA reads to this repeat.
toignore = ["ENSMUSG00000035202"]
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

ISG = utils.ISG.ids()
for file in os.listdir(utils.paths.DE_GENES):
    file = utils.paths.DE_GENES.joinpath(file)
    df = pd.read_csv(file, index_col=0)

    df.index.rename("Ensembl ID", inplace=True)
    names = []
    for geneid in df.index:
        try:
            names.append(utils.ensembl.id_to_name(geneid))
        except ValueError:
            names.append(None)
    df['Gene name'] = names
    df['ISG'] = [x in ISG for x in df.index]
    df.to_csv(file)
