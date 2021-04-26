from multiprocessing import cpu_count


rule index_bam:
    input:
        "{path}.bam"
    output:
        "{path}.bam.bai"
    threads: cpu_count()
    shell:
        "samtools index -@ {threads} {input}"


rule index_fasta:
    input:
        "{path}.fa"
    output:
        "{path}.fa.fai"
    shell:
        "samtools faidx {input}"
