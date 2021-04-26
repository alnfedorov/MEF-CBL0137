from multiprocessing import cpu_count


rule download_sra_pe:
    output:
        r1="results/samples/pe-SRR{ind}/fq.gz/SRR{ind}_1.fq.gz",
        r2="results/samples/pe-SRR{ind}/fq.gz/SRR{ind}_2.fq.gz"
    threads: cpu_count()
    shell:
        "fasterq-dump -o results/samples/pe-SRR{wildcards.ind}/fq.gz/SRR{wildcards.ind}.fq "
        "--progress --threads {threads} SRR{wildcards.ind} && "
        "gzip results/samples/pe-SRR{wildcards.ind}/fq.gz/SRR{wildcards.ind}_1.fq && "
        "gzip results/samples/pe-SRR{wildcards.ind}/fq.gz/SRR{wildcards.ind}_2.fq"
