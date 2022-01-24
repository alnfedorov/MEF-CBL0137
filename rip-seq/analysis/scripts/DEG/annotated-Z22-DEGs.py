import pandas as pd

import utils

merged = None
for ifnb in "48h", "24h", "0h":
    file = utils.paths.DE_GENES.joinpath(f"{ifnb}_Z22_vs_{ifnb}_input.csv")
    df = pd.read_csv(file)
    df = df[(df.padj < utils.DEG.RIP_PADJ) & (df.log2FoldChange > utils.DEG.RIP_LFC)]
    df = df.set_index(['Ensembl Gene ID', 'Gene name', 'Biotype'])[['log2FoldChange']]
    df = df.rename(columns={"log2FoldChange": f"{ifnb}-log2FoldChange"})
    if merged is not None:
        merged = merged.join(df, how='outer')
    else:
        merged = df

df = merged.fillna(0)
df['Max log2FoldChange'] = df.max(axis=1)
assert df['Max log2FoldChange'].isna().sum() == 0
df = df[['Max log2FoldChange']].reset_index()
assert len(df) == len(utils.DEG.Z22())

# Remove IgG enriched genes
df = df[~df['Ensembl Gene ID'].isin(utils.DEG.IgG())]

# ISGs information
ISG = utils.DEG.ISG()
df['ISG'] = df['Ensembl Gene ID'].isin(ISG)

assert ((df['Biotype'] == 'protein_coding') & df['ISG']).sum() == 355

# Editing information
df['Transcript edited in WT (IFNb-72h)'] = df['Gene name'].isin(utils.features.edited_genes())
assert ((df['Biotype'] == 'protein_coding') & df['ISG'] & df['Transcript edited in WT (IFNb-72h)']).sum() == 14

# Transcript features
repeats = pd.read_csv(utils.paths.REPEATS_IN_GENES)
for cls in ['SINE', 'Simple_repeat', 'LTR', 'LINE']:
    subdf = repeats[repeats['Class'] == cls]
    subdf = subdf.groupby(['Gene ID'])['Total repeats'].sum().to_dict()
    df[f"{cls}?"] = df['Ensembl Gene ID'].apply(lambda x: subdf.get(x, 0) > 0)

for group, reps in ["Inverted-SINE",
                    ('Inverted-B4', 'Inverted-Alu', 'Inverted-ID', 'Inverted-B2', 'Inverted-Deu', 'Inverted-MIR')], \
                   ["GU-type", ("(TG)n", "(CA)n")]:
    subdf = repeats[repeats['Name'].isin(reps)]
    subdf = subdf.groupby(['Gene ID'])['Total repeats'].sum().to_dict()
    df[f"{group}?"] = df['Ensembl Gene ID'].apply(lambda x: subdf.get(x, 0) > 0)
df.to_csv(utils.paths.RESULTS.joinpath("Z22-DEGs-details.csv"), index=False)
