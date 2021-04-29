import os
import subprocess
import tempfile

import pandas as pd

from .paths import RESOURCES_FOLDER

ZHUNT_SOURCES = RESOURCES_FOLDER.joinpath("zhunt2.c")
ZHUNT = RESOURCES_FOLDER.joinpath("zhunt2")


def run(query: str, windowsize: int = 8, minsize: int = 6, maxsize: int = 8):
    if not os.path.exists(ZHUNT):
        raise FileNotFoundError(
            f"Compiled zhunt2 is required to predict z-scores. Run the following command to build it: "
            f"gcc {ZHUNT_SOURCES} -lm -o {ZHUNT}")
    fd, temp = tempfile.mkstemp()
    os.close(fd)
    with open(temp, 'w') as stream:
        stream.write(query)

    subprocess.run([ZHUNT, str(windowsize), str(minsize), str(maxsize), temp],
                   check=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, input=query, encoding='ascii')
    with open(temp + ".Z-SCORE", 'r') as stream:
        df = pd.read_csv(stream,
                         names=['Start', 'End', 'not-used-1', 'not-used-2', 'not-used-3',
                                'Z-Score', 'Sequence', 'Conformation'],
                         skiprows=1, sep='\s+')
    os.remove(temp)
    os.remove(temp + ".Z-SCORE")
    return df
