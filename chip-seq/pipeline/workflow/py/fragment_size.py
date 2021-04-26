from typing import List


# fragment size is unlikely to be < 140 for a typical chip-seq
def _parse_macs2_fragment_size(file: str):
    with open(file, 'r') as file:
        results = file.readlines()
        # Search in the output for the # predicted fragment length is X bps
        fragment = [line.strip() for line in results if "# predicted fragment length is" in line]

    if len(fragment) != 1:
        assert "Can't find enough pairs of symmetric peaks to build model!" in results[-2]
        print("Can't determine fragment size, default value is used: 140")
        return 140

    assert len(fragment) == 1, "macs2 predictd changed its output format"
    fragment = fragment[0].split(" ")
    assert fragment[-1] == "bps"

    fragment = int(fragment[-2])
    assert fragment > 0
    return max(140, fragment)


def mean_fragment_size(treatment: List[str]):
    fs = [_parse_macs2_fragment_size(t) for t in treatment]
    return int(round(sum(fs) / len(fs)))
