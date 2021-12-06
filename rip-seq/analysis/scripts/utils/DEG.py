from typing import Set, Optional

import pandas as pd

from . import paths, ensembl

RIP_PADJ = 0.01
RIP_LFC = 1

ISG_PADJ = 0.01
ISG_LFC = 0.585


def _load(file: str, padj: float, lfc: float, biotype: Optional[str] = None) -> Set[str]:
    file = paths.DE_GENES.joinpath(file)
    df = pd.read_csv(file)
    df = df[(df.padj < padj) & (df.log2FoldChange > lfc)]
    if biotype:
        df = df[df.Biotype == biotype]
    return set(df['Ensembl Gene ID'])


def ISG_at_6h(biotype: Optional[str] = None) -> Set[str]:
    return _load("GSE118926.csv", ISG_PADJ, ISG_LFC, biotype)


def ISG_at_24h(biotype: Optional[str] = None) -> Set[str]:
    return _load("IFNb_24h_vs_0h.csv", ISG_PADJ, ISG_LFC, biotype)


def ISG_at_48h(biotype: Optional[str] = None) -> Set[str]:
    return _load("IFNb_48h_vs_0h.csv", ISG_PADJ, ISG_LFC, biotype)


def ISG(biotype: Optional[str] = None) -> Set[str]:
    return set.union(ISG_at_6h(biotype), ISG_at_24h(biotype), ISG_at_48h(biotype))


def IgG_at_0h(biotype: Optional[str] = None) -> Set[str]:
    return _load("0h_IgG_vs_0h_input.csv", RIP_PADJ, RIP_LFC, biotype)


def IgG_at_24h(biotype: Optional[str] = None) -> Set[str]:
    return _load("24h_IgG_vs_24h_input.csv", RIP_PADJ, RIP_LFC, biotype)


def IgG_at_48h(biotype: Optional[str] = None) -> Set[str]:
    return _load("48h_IgG_vs_48h_input.csv", RIP_PADJ, RIP_LFC, biotype)


def IgG(biotype: Optional[str] = None) -> Set[str]:
    return set.union(IgG_at_0h(biotype), IgG_at_24h(biotype), IgG_at_48h(biotype))


def Z22_at_0h(biotype: Optional[str] = None) -> Set[str]:
    return _load("0h_Z22_vs_0h_input.csv", RIP_PADJ, RIP_LFC, biotype)


def Z22_at_24h(biotype: Optional[str] = None) -> Set[str]:
    return _load("24h_Z22_vs_24h_input.csv", RIP_PADJ, RIP_LFC, biotype)


def Z22_at_48h(biotype: Optional[str] = None) -> Set[str]:
    return _load("48h_Z22_vs_48h_input.csv", RIP_PADJ, RIP_LFC, biotype)


def Z22(biotype: Optional[str] = None) -> Set[str]:
    return set.union(Z22_at_0h(biotype), Z22_at_24h(biotype), Z22_at_48h(biotype))


def names(ids: Set[str]) -> Set[str]:
    names = set()
    for ensid in ids:
        try:
            names.add(ensembl.id_to_name(ensid))
        except ValueError:
            continue
    return names
