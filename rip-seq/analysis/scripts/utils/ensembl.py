from .paths import RESOURCES

from functools import lru_cache

import pandas as pd

BIOMART_MAPPING = RESOURCES.joinpath("ENSEMBL-gene-transcripts-info.tsv")


@lru_cache(maxsize=1)
def _load_mapping():
    df = pd.read_csv(BIOMART_MAPPING, sep="\t")
    gene_id_to_gene_name = {k: v for k, v in zip(df['Gene stable ID'].values, df["Gene name"].values)}
    gene_id_to_gene_name.update({k: v for k, v in zip(df['Gene stable ID version'].values, df["Gene name"].values)})
    return gene_id_to_gene_name


def id_to_name(id: str) -> str:
    return _load_mapping()[id]
