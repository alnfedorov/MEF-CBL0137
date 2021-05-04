import os
from pathlib import Path
from subprocess import check_call

import utils

import pandas as pd


os.makedirs(utils.paths.DE_GENES, exist_ok=True)
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
