import pathlib

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import patches as mpatches

import utils

df = pd.read_csv(utils.paths.ADAR1_WT_IFNb_72h_Z22_vs_IgG)
df = df[df['padj'] < 0.1]
assert df['log2FoldChange'].min() >= 0.585
Z22genes = set(df['Ensembl gene ID'])

ISG = set()
for file in [utils.paths.INTERFEROME_MEF, utils.paths.INTERFEROME_ALL]:
    ISG.update(pd.read_csv(file)['Ensembl ID'])
ISG.update(pd.read_excel(utils.paths.ISG_IN_HOME_MICROARRAYS)['ENSEMBLID'])

Z22ISG = Z22genes.intersection(ISG)

keys = ['IFN-I stimulated', "Others"]
shares = [len(Z22ISG), len(Z22genes) - len(Z22ISG)]

colorscheme = {
    "IFN-I stimulated": "#FC766A",
    "Others": "#90BFF9",
}

fig = plt.figure(figsize=(12, 8))
ax = fig.gca()

colors = [colorscheme[k] for k in keys]
wedges, _ = ax.pie(shares, startangle=45,
                   wedgeprops=dict(width=0.5, edgecolor='#585A5B', linewidth=2), colors=colors)

legend = [mpatches.Patch(facecolor=colorscheme[k], label=k, edgecolor='black', linewidth=0.75) for k in keys]
ax.set_xlim(-1, 2.75)
ax.legend(handles=legend, frameon=False, fontsize='xx-large', loc='right')
ax.set_title("Genes bound by Z22 Ab\n(ADAR1-WT, IFNb 72h)", fontsize='xx-large', x=0, y=1.15, transform=ax.transData)
fig.tight_layout()

saveto = pathlib.Path(__file__).name.replace(".py", ".eps")
fig.savefig(utils.paths.RESULTS.joinpath(saveto), bbox_inches="tight", pad_inches=0)
