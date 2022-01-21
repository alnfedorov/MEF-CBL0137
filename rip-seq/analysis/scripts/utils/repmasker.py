import os
from functools import lru_cache
from typing import Set

from .paths import RESOURCES


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


_REPEAT_MASKER_CLASSIFICATION = os.path.join(RESOURCES, 'repeat-masker-name-class-family.tsv')


def classify(repname: str):
    classification = _load_classification(_REPEAT_MASKER_CLASSIFICATION)
    return classification[repname]


@lru_cache(maxsize=1)
def family_to_cls(family: str) -> str:
    mapping = {x[1]: x[2] for x in _load_classification(_REPEAT_MASKER_CLASSIFICATION).values()}
    return mapping[family]


@lru_cache(maxsize=1)
def names() -> Set[str]:
    return set(x[0] for x in _load_classification(_REPEAT_MASKER_CLASSIFICATION).values())


@lru_cache(maxsize=1)
def families() -> Set[str]:
    return set(x[1] for x in _load_classification(_REPEAT_MASKER_CLASSIFICATION).values())


@lru_cache(maxsize=1)
def classes() -> Set[str]:
    return set(x[2] for x in _load_classification(_REPEAT_MASKER_CLASSIFICATION).values())
