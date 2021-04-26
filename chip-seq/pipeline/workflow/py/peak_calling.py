class peak_calling:
    class macs2:
        @staticmethod
        def flags(design: str) -> str:
            params = config['peak_calling']['macs2']

            flags = []
            if "gsize" in params:
                flags.append(f"-g {params['gsize']}")
            if "fdr" in params:
                flags.append(f"-q {params['fdr']}")
            if "fecutoff" in params:
                flags.append(f"--fe-cutoff={params['fecutoff']}")
            if "keepdup" in params:
                flags.append(f"--keep-dup={params['keepdup']}")

            if Design(design).ispe:
                flags.append(" --format=BAMPE")
            else:
                flags.append("--nomodel")
                flags.append(f"--extsize={mean_fragment_size(input.fragment_sizes)}")

            return flags
