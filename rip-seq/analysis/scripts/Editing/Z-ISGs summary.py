import pandas as pd

import utils

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

REP_CLS_TO_REPORT = {"SINE", "Low_complexity", "Simple_repeat", "LTR", "LINE"}

ZRNA = utils.DEG.Z22('protein_coding') - utils.DEG.IgG()
ISG = utils.DEG.ISG()

# Double check that all edits in coding Z-ISGs happen in 3`UTRs
df = pd.read_csv(utils.paths.RESULTS.joinpath("editing-distribution.csv"))
ZISGs = utils.DEG.names(ISG & ZRNA)
df = df[df['Gene'].apply(lambda x: any(gene in ZISGs for gene in x.split(';')))]
assert (df['Location'] == 'utr3').all()

# Get all edited repeats
EDITED = utils.features.edited_repeats(utils.DEG.names(ZRNA))

# Download a summary table
df = pd.read_csv(utils.paths.REPEATS_IN_GENES)
df = df[df['Gene ID'].isin(ZRNA)]
assert df['Total repeats'][df['Edited in ADAR1 WT']].sum() == len(EDITED)
df['ISG'] = df['Gene ID'].isin(ISG)

groupas = ['Class', 'Name', 'ISG', 'Edited in ADAR1 WT']
# Report selected records separately
separate = {"(TG)n", "(CA)n",
            'Inverted-Alu', 'Inverted-B2', 'Inverted-B4', 'Inverted-Deu', 'Inverted-ID', 'Inverted-MIR'}

# Create a gene-level summary
genes = df.groupby(['Gene ID', 'ISG', 'Class', 'Edited in ADAR1 WT']).sum()
assert (genes['Total repeats'] > 0).all()
genes['Name'] = 'Total'
genes = genes.reset_index().drop(columns=['Gene ID']).groupby(groupas).count()

selgenes = df[df['Name'].isin(separate)]
assert (selgenes['Total repeats'] > 0).all()
selgenes = selgenes.drop(columns=['Gene ID']).groupby(groupas).count()

genes = pd.concat([genes, selgenes]).rename(columns={"Total repeats": "Total genes"}).reset_index()
genes = genes[genes['Class'].isin(REP_CLS_TO_REPORT)].set_index(groupas).sort_index()
print("Summary of coding Z-genes with a given feature")
print(genes)

# Create a repeat-level summary
df = df.drop(columns=['Gene ID'])
df = df.groupby(groupas).sum().reset_index()
separate = df[df['Name'].isin(separate)]

df["Name"] = "Total"
df = df.groupby(groupas).sum()
separate = separate.groupby(groupas).sum()

df = pd.concat([df, separate]).reset_index()
assert df[df['Edited in ADAR1 WT'] & (df["Name"] == "Total")]['Total repeats'].sum() == len(EDITED)

# Print results for selected classes
df = df[df['Class'].isin(REP_CLS_TO_REPORT)]
df = df.set_index(groupas).sort_index()
print("Number of repeats in each class located in coding Z-RNAs")
print(df)
