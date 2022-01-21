from pybedtools import BedTool

import utils

repgroups = [
    "LINE", "Low_complexity", "LTR", "RC", "RNA", "Satellite", "scRNA",
    "Simple_repeat", "snRNA", "srpRNA"
]

repgroups = {k: [] for k in repgroups}

for repeat in BedTool(utils.paths.REPMASKER):
    _, _, group = utils.repmasker.classify(repeat.name)
    if group in repgroups:
        repgroups[group].append(repeat)

utils.paths.EDITING_INDEX_INTERVALS.mkdir(parents=True, exist_ok=True)
for name, intervals in repgroups.items():
    bed = BedTool(intervals).sort()
    bed.saveas(
        utils.paths.EDITING_INDEX_INTERVALS.joinpath(f"{name}.bed").absolute()
    )
