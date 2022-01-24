import scipy.stats

from utils import features, DEG

background = features.expressed_genes(gtype='protein_coding')
ZRNA = DEG.Z22(biotype='protein_coding') - DEG.IgG()
ISG = DEG.ISG(biotype='protein_coding')
assert ZRNA.issubset(background) and ISG.issubset(background)

# Universe = background
# Feature = ZRNA
# Sample = ISG

bckgall = len(background)  # Total size of the background
bckgev = len(ZRNA)  # Number of records with a given feature

smplall = len(ISG)  # Total size of the sample
smplev = len(ISG & ZRNA)  # Number of records in the sample with a given feature
smplnonev = smplall - smplev

# [`M`, `n`, `N`]
# M - total number of objects
# n - total number of objects with a given feature
# N - sample size
pv = scipy.stats.hypergeom(bckgall, bckgev, smplnonev + smplev).sf(smplev - 1)

print(f"Total expressed coding genes: {len(background)}")
print(f"Total coding Z-RNAs: {len(ZRNA)}")
print(f"Total coding ISGs: {len(ISG)}")
print(f"Coding ISGs & Z-RNA: {len(ISG & ZRNA)}")
print(f"Hypergeometric mean p-value: {pv}")
