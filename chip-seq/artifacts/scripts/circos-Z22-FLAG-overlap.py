import os
import pathlib
from pathlib import Path
from pybedtools import BedTool
from subprocess import check_call

import utils

blacklisted = BedTool(utils.paths.BLACKLIST)

z22 = BedTool(utils.paths.Z22_CURAX_14H).subtract(blacklisted, A=True)
flag = BedTool(utils.paths.FLAG_CURAX_14H).subtract(blacklisted, A=True)
intersection = z22.intersect(flag)

z22, flag, intersection = z22.saveas().fn, flag.saveas().fn, intersection.saveas().fn

script = Path(__file__).absolute().as_posix().replace(".py", ".R")

saveto = pathlib.Path(__file__).name.replace(".py", ".eps")
saveto = utils.paths.RESULTS.joinpath(saveto)

check_call(["Rscript", script, z22, flag, intersection, saveto])
os.unlink(z22)
os.unlink(flag)
os.unlink(intersection)
