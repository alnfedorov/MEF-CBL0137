import os
from functools import lru_cache

from .paths import RESOURCES_FOLDER

L1Md_A_CONSENSUS = RESOURCES_FOLDER.joinpath("L1Md-A-consensus.txt")
L1Md_T_CONSENSUS = RESOURCES_FOLDER.joinpath("L1Md-T-consensus.txt")


@lru_cache(2)
def sequence(l1: str) -> str:
    path = {"L1Md_A": L1Md_A_CONSENSUS, "L1Md_T": L1Md_T_CONSENSUS}[l1]
    assert os.path.exists(path), "run build-L1Md-A-T-consensus.py to generate consensus sequences for L1Md_A and L1Md_T"
    with open(path, 'r') as stream:
        seq = stream.read().strip().upper()
    return seq


# created by the NCBI ORFfinder online tool
def structure(l1: str) -> dict:
    return {
        "L1Md_T": {
            "5`UTR": (0, 901),
            "ORF1": (901, 2016),
            "ORF2": (2057, 5902),
            "3`UTR": (5902, len(sequence('L1Md_T')))
        },
        "L1Md_A": {
            "5`UTR": (0, 926),
            "ORF1": (926, 1999),
            "ORF2": (2040, 5885),
            "3`UTR": (5885, len(sequence('L1Md_A')))
        },
    }[l1]
