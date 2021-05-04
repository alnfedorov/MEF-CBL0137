import os
from io import StringIO

import pandas as pd

import utils

# Merge all files into a single string
alldata = [
    "Dataset ID	Fold Change	Inteferome Type	Treatment Time	Gene Name	Description	GenBank Accession	Ensembl ID	Probe ID\n"
]
folder = utils.paths.INTERFEROME_ALL_QUERIES
for query in os.listdir(folder):
    file = os.path.join(folder, query, "results.txt")
    with open(file, 'r') as file:
        skipped = 0
        for line in file:
            if line == alldata[0]:
                break
            skipped += 1
        print("Skipped header lines: ", skipped)
        alldata.extend(file)
alldata = "".join(alldata)

df = pd.read_csv(StringIO(alldata), sep="\t")
df = df[['Fold Change', 'Inteferome Type', 'Gene Name', 'Ensembl ID']]

# Filter by median fold change
df = df.groupby(['Inteferome Type', 'Gene Name', 'Ensembl ID']).median().reset_index()
df = df.rename(columns={"Fold Change": "Median fold change"})
df = df[(df['Median fold change'] > 1.5) & (df['Inteferome Type'] == 'I')]

# Save results
df = df[['Gene Name', 'Ensembl ID', 'Median fold change']]
df.to_csv(utils.paths.INTERFEROME_ALL, index=False)

# similar for the MEF
df = pd.read_csv(utils.paths.INTERFEROME_MEF_QUERY, skiprows=19, sep='\t')
df = df[['Fold Change', 'Inteferome Type', 'Gene Name', 'Ensembl ID']]
df = df.groupby(['Inteferome Type', 'Gene Name', 'Ensembl ID']).median().reset_index()
df = df.rename(columns={"Fold Change": "Median fold change"})
df = df[(df['Median fold change'] > 1.5) & (df['Inteferome Type'] == 'I')]
df.to_csv(utils.paths.INTERFEROME_MEF, index=False)
