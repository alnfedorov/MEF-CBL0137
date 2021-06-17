from pybedtools import BedTool, Interval

import utils

covblacklist = {
    "7SK", "5S", "LSU-rRNA_Hsa", "SSU-rRNA_Hsa", "7SLRNA"
}

bed = BedTool(utils.paths.REPMASKER).filter(lambda x: "tRNA" in x.name or x.name in covblacklist).sort()
bed = bed.slop(b=250, genome='mm10').sort().merge()
bed = list(bed)

# rRNA genes
bed.append(Interval("chr4", 43492707, 43493127))
# bed.append(Interval("chrM", 0, 1036))
bed.append(Interval("chr17", 39842858, 39846364))
bed = BedTool(bed).sort().merge()
bed.saveas(utils.paths.COVERAGE_BLACKLIST)
