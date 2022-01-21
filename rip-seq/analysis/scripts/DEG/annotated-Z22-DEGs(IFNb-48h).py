import pandas as pd

import utils

ISG = utils.DEG.ISG()

file = utils.paths.DE_GENES.joinpath("48h_Z22_vs_48h_input.csv")
df = pd.read_csv(file).drop(columns=['Z22-DEG'])
df['Transcript edited in WT (IFNb-72h)'] = df['Ensembl Gene ID'].isin(utils.features.edited_genes())
assert (df['ISG'] == df['Ensembl Gene ID'].isin(utils.DEG.ISG())).all()

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
df.to_csv(utils.paths.RESULTS.joinpath("Z22-DEG(IFNb-48h)-details.csv.gz"), index=False)
