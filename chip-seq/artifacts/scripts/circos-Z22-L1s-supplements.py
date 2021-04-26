import os
import pathlib
from pathlib import Path
from subprocess import check_call

from pybedtools import BedTool

import utils

blacklisted = BedTool(utils.paths.BLACKLIST)
z22 = BedTool(utils.paths.Z22_CURAX_14H).subtract(blacklisted, A=True)

repmasker = BedTool(utils.paths.REPEATMASKER)
repmasker = repmasker.subtract(blacklisted, A=True).intersect(z22, u=True, wa=True)

l1s = {"L1Md_A": [], "L1Md_T": []}
for repeat in repmasker:
    if repeat.name in l1s:
        l1s[repeat.name].append(repeat)
l1s = {k: BedTool(v).sort().saveas().fn for k, v in l1s.items()}

l1mda, l1mdt = l1s['L1Md_A'], l1s['L1Md_T']

script = Path(__file__).absolute().as_posix().replace(".py", ".R")

saveto = pathlib.Path(__file__).name.replace(".py", ".eps")
saveto = utils.paths.RESULTS.joinpath(saveto)

check_call(["Rscript", script, l1mdt, l1mda, saveto])
os.unlink(l1mdt)
os.unlink(l1mda)
