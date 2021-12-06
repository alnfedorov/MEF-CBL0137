from pathlib import Path
from subprocess import check_call

import pandas as pd

import utils

utils.paths.DE_GENES.mkdir(parents=True, exist_ok=True)

# GSE118926
script = Path(__file__).absolute().as_posix().replace(".py", "(GSE118926).R")
folder, saveto = utils.paths.SALMON.joinpath("GSE118926"), utils.paths.DE_GENES
check_call(["Rscript", script, folder, saveto])
GSE118926_degs = utils.paths.DE_GENES.joinpath("GSE118926.alldegs.csv")
assert utils.paths.DE_GENES.joinpath("GSE118926.csv").is_file() and \
       GSE118926_degs.is_file()

# Our data
script = Path(__file__).absolute().as_posix().replace(".py", "(ADAR1-KO_ZBP1-KO).R")
folder, saveto = utils.paths.SALMON.joinpath("ADAR1-KO_ZBP1-KO"), utils.paths.DE_GENES
check_call(["Rscript", script, folder, saveto, GSE118926_degs])

# Postprocess and annotate
for file in utils.paths.DE_GENES.glob("*.csv"):
    if file.name in ("expression.normalized.csv", "IFNb: cluster.csv"):
        continue

    df = pd.read_csv(file, index_col=0)

    names, biotypes = [], []
    df.index.rename("Ensembl Gene ID", inplace=True)
    for geneid in df.index:
        try:
            names.append(utils.ensembl.id_to_name(geneid))
            biotypes.append(utils.ensembl.gene_id_to_type(geneid))
        except ValueError:
            names.append(None)
            biotypes.append(None)
    df['Gene name'] = names
    df['Biotype'] = biotypes

    file.unlink()
    file = file.as_posix().replace("all_", "").replace("selection_", "")
    df.to_csv(file)

ISG = utils.DEG.ISG()
Z22 = utils.DEG.Z22()
for file in utils.paths.DE_GENES.glob("*.csv"):
    if file.name in ("expression.normalized.csv", "IFNb: cluster.csv"):
        continue
    df = pd.read_csv(file)

    df['ISG'] = df['Ensembl Gene ID'].isin(ISG)
    df['Z22-DEG'] = df['Ensembl Gene ID'].isin(Z22)
    df.to_csv(file, index=False)
