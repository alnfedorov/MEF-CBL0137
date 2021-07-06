import pandas as pd
import scipy.stats as stats

import utils

# 3’UTRs of ISG mRNAs were disproportionately targeted for editing by ADAR1 p150 (hypergeometric)
# Backrground - all annotated 3`UTRs
df = pd.read_csv(
    utils.paths.RESULTS.joinpath("RAW-Editing-distribution(Z22 & editing index > 0.5% & mean coverage > 10).csv")
)
df = df[(df['Location'] == "utr3")]

edited_allutr3 = len(df)  # all tries
isgutr3 = len(df[df['ISG']])  # successes

isgs = utils.ISG.names()

# fetch all annotated 3`UTRs
set(utils.ensembl.transcript_to_gene_name(interval.name.split("_")[0]) for interval in utr3)

