import shutil
from functools import reduce
from pathlib import Path
from subprocess import run

import pandas as pd

import utils

folder = Path("/tmp/MEF-ISG")
folder.mkdir(exist_ok=True)

# download and unpack data
run(f'tar xv -C {folder} < <(wget -q -O - "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE128110&format=file")',
    check=True, shell=True, executable='/bin/bash')
run(f"cd {folder} && gunzip * && rm -rf *IRF9ko*", check=True, shell=True)

# Build meta and counts tables
title, filename, condition = [], [], []
for file in folder.iterdir():
    if file.name == "samples.tsv" or file.name == "counts.tsv":
        continue
    filename.append(file.name)
    title.append("_".join(file.name.split("_")[1:3]))
    title[-1] = title[-1].replace("-", "_")
    if "_b-" in file.name:
        condition.append("IFNb")
    else:
        assert "_ut-" in file.name
        condition.append("untreated")

samples = pd.DataFrame({"title": title, "filename": filename, "condition": condition}).set_index("title")
samples.to_csv(folder.joinpath("samples.tsv"), sep="\t", index=True)

counts = []
for title, file in zip(samples.index, samples['filename']):
    df = pd.read_csv(folder.joinpath(file), sep="\s+").set_index('gene_id')
    df.rename(columns={'read_counts': title}, inplace=True)
    counts.append(df)

counts = reduce(lambda left, right: left.join(right), counts)
counts.to_csv(folder.joinpath("counts.tsv"), sep="\t", index=True)

# run DEseq2 to infer ISG
script = Path(__file__).absolute().as_posix().replace(".py", ".R")
run(["Rscript", script, folder, utils.paths.ISG_MEF_NGS], check=True)
shutil.rmtree(folder)

df = pd.read_csv(utils.paths.ISG_MEF_NGS, index_col=0)
df = df[df['padj'] < 0.1]
df.index.rename("Gene name", inplace=True)

ensemblid = []
for gname in df.index:
    gid = utils.ensembl.name_to_id(gname)
    ensemblid.append(gid)

df['Ensembl ID'] = ensemblid
df = df.reset_index()[['Gene name', 'Ensembl ID', 'log2FoldChange', 'padj']]
df.to_csv(utils.paths.ISG_MEF_NGS, index=False)
