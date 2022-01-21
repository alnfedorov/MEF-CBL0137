from collections import defaultdict
from itertools import chain

import pandas as pd

import utils
from utils.miscellaneous import hashinter

genes = utils.features.expressed_genes()
assert utils.DEG.Z22().issubset(genes) and utils.DEG.ISG().issubset(genes)

transcripts = {k: utils.ensembl.gene_id_to_transcripts(k) for k in genes}
transcripts = set(chain(*transcripts.values()))

# Extract all relevant exons
EXONS = utils.features.exons(transcripts)

# Use mRNA annotation for coding genes, which is more complete
mrna = utils.features.mRNA(transcripts)
for k, v in mrna.items():
    if k in EXONS:
        assert sum(x.length for x in EXONS[k]) <= sum(x.length for x in v)
    EXONS[k] = v

# Get all edited repeats
EDITED = utils.features.edited_repeats(utils.DEG.names(genes))

# Build a summary table
repeats = utils.features.repeats_in(EXONS)
data = defaultdict(int)  # Total number of repeats hosted by a gene
counted = defaultdict(
    set)  # Track counted repeats per gene and don't count them twice (overlap between multiple isoforms)

# Count long Z-prone isoforms first
ZRNA = [utils.ensembl.gene_id_to_transcripts(x) for x in utils.DEG.Z22() - utils.DEG.IgG()]
ZRNA = set(chain(*ZRNA))
order = {k: (k in ZRNA, sum(x.length for x in v)) for k, v in EXONS.items()}
order = sorted(order.items(), key=lambda x: x[1], reverse=True)

for tid, _ in order:
    if tid not in repeats:
        continue

    reps = repeats[tid]
    gid = utils.ensembl.transcript_to_gene_id(tid)

    sines = defaultdict(list)  # match inverted SINEs from the same family later
    for r in reps:
        _, family, cls = utils.repmasker.classify(r.name)
        family, cls = family.replace("?", ""), cls.replace("?", "")

        hash = hashinter(r, cls)
        if hash in counted[gid]:
            continue

        if cls == "SINE":
            sines[family].append(r)

        counted[gid].add(hash)
        data[(gid, r.name, cls, hash in EDITED)] += 1

    for family, reps in sines.items():
        matched, singletons = utils.features.match_inverted(reps, EDITED)
        for s1, s2 in matched:
            for x in s1, s2:
                hash = hashinter(x, "SINE")
                edited = hash in EDITED
                assert (gid, x.name, "SINE", edited) in data
                data[(gid, x.name, "SINE", edited)] -= 1
                assert data[(gid, x.name, "SINE", edited)] >= 0
                data[(gid, f"Inverted-{family}", "SINE", edited)] += 1

# All edited repeats must be counted
counted = set(chain(*counted.values()))
assert all(x in counted for x in EDITED), EDITED - counted

# Build a final DataFrame
key, values = zip(*data.items())
gid, name, cls, isedited = zip(*key)
df = pd.DataFrame(
    {"Gene ID": gid, "Class": cls, "Name": name, "Edited in ADAR1 WT": isedited, "Total repeats": values}
)
# One edited repeat can be shared between multiple genes
assert df['Total repeats'][df['Edited in ADAR1 WT']].sum() >= len(EDITED)
assert (df['Total repeats'] >= 0).all()
df.to_csv(utils.paths.REPEATS_IN_GENES, index=False)
