import math

import pybedtools

import utils

Q_VALUES_COLUMN = 8
FC_COLUMN = 6

blacklisted = pybedtools.BedTool(utils.paths.BLACKLIST).sort()

for name, bed in [("0h: FLAG vs IgG", utils.paths.FLAG_CURAX_0H),
                  ("14h: FLAG vs IgG", utils.paths.FLAG_CURAX_14H),
                  ("14h: Z22 vs IgG", utils.paths.Z22_CURAX_14H)]:
    print(name)
    bed = pybedtools.BedTool(bed).subtract(blacklisted, A=True)

    fdr005 = len([x for x in bed
                  if float(x.fields[Q_VALUES_COLUMN]) > -math.log10(0.05) and float(x.fields[FC_COLUMN]) > 2])
    assert fdr005 == len(bed)
    print(f"\tFDR<0.05 & FC>2: {fdr005}")

    fe5 = len([x for x in bed if float(x.fields[FC_COLUMN]) > 5])
    print(f"\tFDR<0.05 & FC>5: {fe5}")
