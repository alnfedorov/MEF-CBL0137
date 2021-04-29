import os
from functools import lru_cache

from .paths import RESOURCES_FOLDER


@lru_cache(maxsize=1)
def _load_classification(path: str):
    mapping = {}
    with open(path, 'r') as stream:
        for line in stream:
            line = line.strip()
            if len(line) == 0:
                continue
            name, cls, family = line.split()
            assert name not in mapping, name
            mapping[name] = (name, family, cls)
    return mapping


_REPEAT_MASKER_CLASSIFICATION = os.path.join(RESOURCES_FOLDER, 'repeat-masker-name-class-family.tsv')


def classify(repname: str):
    classification = _load_classification(_REPEAT_MASKER_CLASSIFICATION)
    return classification[repname]
