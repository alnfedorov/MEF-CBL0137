import subprocess
from functools import reduce
from multiprocessing import cpu_count
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd

import utils


def findfile(folder: Path, prefix: str, postfix: Optional[str] = None) -> Path:
    files = []
    for file in folder.iterdir():
        if file.name.startswith(prefix) and (postfix is None or file.name.endswith(postfix)):
            files.append(file)
    assert len(files) == 1, f"{prefix} -> {files}"
    return files[0]


def calculate_scaling(samples: List[str]) -> List[float]:
    """Median of ratios approach"""
    files = [findfile(utils.paths.SALMON, sample) for sample in samples]
    assert all(x.is_file() for x in files)

    frames = []
    for key, file in zip(samples, files):
        df = pd.read_csv(file.absolute(), sep="\t")[['Name', 'NumReads']].set_index('Name')
        frames.append(
            df.rename(columns={"NumReads": key})
        )
    df = reduce(lambda left, right: left.join(right), frames)
    df = df[(df > 0).all(axis=1)]
    df['reference'] = np.exp(df.apply(np.log).mean(axis=1))

    scaling = []
    for key in samples:
        ratios = df[key] / df['reference']
        scaling.append(1 / float(np.median(ratios)))
    return scaling


def bamCompare(treatment: Path, treatment_scale: float, control: Path, control_scale: float, saveto: str):
    subprocess.run([
        "bamCompare", "-b1", treatment, "-b2", control, "-o", saveto, "--operation", "ratio",
        "--scaleFactors", f"{treatment_scale}:{control_scale}", "--verbose", "-p", str(cpu_count() - 2),
        "--samFlagExclude", "2816"  # not primary + failed QC check + supplementary
    ], check=True)


allsamples = [x.name.rstrip(".quant.genes.sf") for x in utils.paths.SALMON.iterdir()]
for cell in ("ADAR1-WT", "ADAR1-KO"):
    for ifnb in ("IFNb-0h", "IFNb-72h"):
        if cell == "ADAR1-KO":
            antibodies = ["IgG", "J2", "Z22", "FLAG"]
        else:
            antibodies = ["IgG", "J2", "Z22"]

        keys = [f"{cell}-{ifnb}-{ab}" for ab in antibodies]
        samples = [sample for sample in allsamples if any(k in sample for k in keys)]

        # normalized counts = scale * counts
        scaling = calculate_scaling(samples)
        scaling = {k: v for k, v in zip(samples, scaling)}

        # experiment id
        sortkey = lambda x: x.split("-")[1][1:]
        controls = sorted([x for x in samples if "IgG" in x], key=sortkey)
        for ab in antibodies:
            if ab == "IgG":
                continue
            treatments = sorted([x for x in samples if ab in x], key=sortkey)
            for replicate, (treatment, control) in enumerate(zip(treatments, controls)):
                treatmentbam = findfile(utils.paths.BAM, treatment, ".bam")
                controlbam = findfile(utils.paths.BAM, control, ".bam")

                saveto = utils.paths.SIGNAL.joinpath(f"{cell}-{ifnb}-{ab}-vs-IgG-{replicate}.bw")
                bamCompare(treatmentbam, scaling[treatment], controlbam, scaling[control], saveto)

        # for sample, scale in zip(samples, scaling):
        #     bam = findfile(utils.paths.BAM, sample, ".bam")
        #     saveto = utils.paths.SIGNAL.joinpath(sample + ".bw")
        #     subprocess.check_call([
        #         "bamCoverage", "-p", str(cpu_count()), "-b", str(bam.absolute()), "-o", str(saveto.absolute()),
        #         f"--scaleFactor", str(scale), "--verbose"
        #     ])
