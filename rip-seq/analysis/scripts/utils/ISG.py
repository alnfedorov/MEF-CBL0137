from typing import Set

import pandas as pd

from . import paths, ensembl


def ids(padjthr: float = 0.1) -> Set[str]:
    ISG = set()

    # Interferome median > 1.5
    ISG.update(pd.read_csv(paths.ISG_INTERFEROME)['Ensembl ID'])

    # MEF specific
    ISG.update(pd.read_excel(paths.ISG_MEF_MICROARRAY)['ENSEMBLID'])

    df = pd.read_csv(paths.ISG_MEF_NGS)
    df = df[df['padj'] < padjthr]
    assert df['log2FoldChange'].min() > 0.585
    ISG.update(df['Ensembl ID'])
    return ISG


def names(padjthr: float = 0.1) -> Set[str]:
    names = set()
    for ensid in ids(padjthr):
        try:
            names.add(ensembl.id_to_name(ensid))
        except ValueError:
            continue
    return names
