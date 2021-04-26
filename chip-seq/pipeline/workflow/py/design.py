from typing import List


class Design:
    def __init__(self, design):
        self.title = design
        design = config['design'][design]

        self.treatment = design["treatment"]
        self.control = design["control"]

        self.ispe = samples.ispaired(self.treatment[0])
        assert all(samples.ispaired(x) == self.ispe for x in self.treatment + self.control)

    def _files(self, allsamples: List[str], mode: str):
        return [samples.bam(s, mode) for s in allsamples]

    def dependencies(self, mode: str):
        inputs = {
            "treatment": self.treatment_files(mode), 
            "control": self.control_files(mode),
            "indices": self.index_files(mode)
        }
        if not self.ispe:
            inputs["fragment_sizes"] = self.fragment_size_files()
        return inputs

    def treatment_files(self, mode: str):
        return self._files(self.treatment, mode)

    def control_files(self, mode: str):
        return self._files(self.control, mode)

    def index_files(self, mode: str):
        bam = self.treatment_files(mode) + self.control_files(mode)
        return [b + ".bai" for b in bam]

    def fragment_size_files(self):
        if self.ispe:
            raise ValueError("fragment size estimation is not necessary for the PE fragments")
        return [f"results/samples/se-{sample}/qc/macs2-fragment-size.txt" for sample in self.treatment]
