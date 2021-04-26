from typing import List


# Namespace-like
class samples:
    @staticmethod
    def all() -> List[str]:
        return config["samples"]["paired"]
 
    @staticmethod
    def ispaired(sample: str) -> bool:
        return sample in config["samples"]["paired"]

    @staticmethod
    def folder(sample: str) -> str:
        prefix = "pe" if samples.ispaired(sample) else "se"
        return f"results/samples/{prefix}-{sample}"

    @staticmethod
    def bam(sample: str, mode: str = "clean") -> str:
        assert mode in ("clean", "raw")        
        postfix = {
            "clean": "-clean.sorted.dupsmarked.bam",
            "raw": ".sorted.dupsmarked.bam"
        }[mode]
        return f"{samples.folder(sample)}/bam/{sample}{postfix}"

    @staticmethod
    def bamindex(sample: str, mode: str = "clean") -> str:
        return samples.bam(sample, mode) + ".bai"

    @staticmethod
    def organism(sample: str) -> str:
        for organism, v in config["samples"]["organism"].items():
            if sample in v:
                return organism
