import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import utils

plt.rcParams['svg.fonttype'] = 'none'

# Download mapping, gene -> cluster / z-scores
clusters = pd.read_csv(utils.paths.DE_GENES.joinpath("IFNb: clusters.csv"), index_col=0)
clusters = clusters[['genes', 'value', 'cutoff0.444', 'IFNb']]
clusters.rename(columns={"value": "Z-score", "IFNb": "IFNβ", "cutoff0.444": "cluster"}, inplace=True)

# Print target genes cluster and save per-cluster genes
symbols = {k: utils.DEG.names(v) for k, v in clusters.groupby('cluster').genes}
for gene in ["Knl1", "Eif2ak2", "Ifih1", "Slfn5", "Ddx58", "Tapbp", "Agl", "Xrn1"]:
    ids = [k for k, v in symbols.items() if gene in v]
    print(f"{gene} -> {ids}")

utils.paths.DE_CLUSTERS.mkdir(parents=True, exist_ok=True)
for clind, df in clusters.groupby('cluster'):
    with open(utils.paths.DE_CLUSTERS.joinpath(f"cluster-{clind}.txt"), 'w') as stream:
        for gene in df.genes.unique():
            print(gene, file=stream)

# Get the estimated expression at each timepoint
expression = pd.read_csv(utils.paths.DE_GENES.joinpath("expression.normalized.csv"), index_col=0).T.reset_index()
expression['IFNβ'] = expression['index'].apply(lambda x: x.split('_')[2].lstrip("IFNb-"))
expression.drop(columns='index', inplace=True)
expression = expression.groupby('IFNβ').mean().T

# Cumulative expression at each timepoint
TOTAL_EXPRESSION = expression.sum().to_dict()
expression = expression.to_dict()

# Annotated Z-formers
Z22 = utils.DEG.Z22().difference(utils.DEG.IgG())
clusters['Z-former'] = clusters['genes'].isin(Z22)

# Correct timepoints order
order = ["0h", "24h", "48h"]
clusters = clusters.sort_values("IFNβ", key=lambda x: x.apply(lambda y: order.index(y)))

assert len(clusters.cluster.unique()) == 6

# Save two versions of the plot, with and without per-gene lines
pallette = sns.palettes.SEABORN_PALETTES['muted']
for lines_only in True, False:
    fig, axes = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(16, 10))
    axes = axes.ravel()

    # iterate clusters
    for clind, ax in enumerate(axes):
        subdf = clusters[clusters.cluster == clind + 1]
        # plot per-gene lines in each cluster if needed
        if lines_only:
            subdf.groupby('genes').apply(
                lambda x: ax.plot(x['IFNβ'], x['Z-score'], color=pallette[clind], linewidth=0.5, alpha=0.05)
            )
        # boxplot
        sns.boxplot(data=subdf, x='IFNβ', y='Z-score', color=pallette[clind], ax=ax, order=order, fliersize=3)

        # Total number of genes
        degs = len(subdf.genes.unique())
        ax.set_title(f"$\\bf{{Cluster~{clind + 1}}}$\n({degs} DEGs)", loc='right')
        # Total number of Z-formers
        zformers = subdf['Z-former'].sum()
        text = f"(cluster {zformers / degs * 100:.2f}% / overall {zformers / len(Z22) * 100:.2f}%)"
        ax.set_title(f"Z22 enriched: {zformers}\n{text}", loc='left', fontsize=11)

        # Print cluster expression only for version without per-gene lines
        if not lines_only:
            for ind, ifnb in enumerate(order):
                # cumulative expression
                expr = subdf.genes.apply(lambda x: expression[ifnb][x]).sum() / TOTAL_EXPRESSION[ifnb]
                # cumulative expression for Z-genes
                zexpr = subdf[subdf['Z-former']].genes.apply(
                    lambda x: expression[ifnb][x]
                ).sum() / TOTAL_EXPRESSION[ifnb]
                # label
                text = f"$\\bf{{∀:}}$ {expr * 100:.2f}%\n$\\bf{{Z:}}$ {zexpr * 100:.2f}%"
                ax.text(
                    0.01, 0.98 - ind * 0.15, text,
                    transform=ax.transAxes, horizontalalignment='left', verticalalignment="top", fontsize=11
                )
    if lines_only:
        fig.savefig(utils.paths.RESULTS.joinpath("IRGs clusters.lines.svg"), bbox_inches='tight', pad_inches=0)
    else:
        fig.savefig(utils.paths.RESULTS.joinpath("IRGs clusters.svg"), bbox_inches='tight', pad_inches=0)
