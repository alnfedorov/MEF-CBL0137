from multiprocessing import cpu_count


rule trim_galore_pe:
    input:
        r1="results/samples/pe-{sample}/fq.gz/{sample}_1.fq.gz",
        r2="results/samples/pe-{sample}/fq.gz/{sample}_2.fq.gz"
    output:
        trimmed_r1="results/samples/pe-{sample}/fq.gz/{sample}_1_val_1.fq.gz",
        trimmed_r2="results/samples/pe-{sample}/fq.gz/{sample}_2_val_2.fq.gz",
        report_r1="results/samples/pe-{sample}/fq.gz/{sample}_1.fq.gz_trimming_report.txt",
        report_r2="results/samples/pe-{sample}/fq.gz/{sample}_2.fq.gz_trimming_report.txt"
    threads: min(cpu_count(), 4)
    shell:
        "trim_galore --stringency 5 --cores {threads} --paired --gzip -o results/samples/pe-{wildcards.sample}/fq.gz/ {input.r1} {input.r2}"


rule minimap2_align_pe:
    input:
        r1="results/samples/pe-{sample}/fq.gz/{sample}_1_val_1.fq.gz",
        r2="results/samples/pe-{sample}/fq.gz/{sample}_2_val_2.fq.gz",
        index="resources/mm10/genome.mmi"
    output:
        bam=temp("results/samples/pe-{sample}/bam/{sample}.bam")
    threads: cpu_count()
    shell: "minimap2 -t {threads} -ax sr {input.index} {input.r1} {input.r2} | samtools view -b > {output.bam}"


rule markdups_pe:
    input:
        "results/samples/pe-{sample}/bam/{sample}.bam"
    output:
        "results/samples/pe-{sample}/bam/{sample}.sorted.dupsmarked.bam"
    threads: cpu_count()
    shell: # picard fails with heap out of memory even for ~100gb JVM heap size
        "samtools sort -@ {threads} -n {input} | "
        "samtools fixmate -@ {threads} -m - - | "
        "samtools sort -@ {threads} | "
        "samtools markdup -S -@ {threads} -s - - | "
        "samtools sort -@ {threads} -l9 -o {output} -"


rule clean_reads_pe:
    input:
        "results/samples/pe-{sample}/bam/{sample}.sorted.dupsmarked.bam"
    output:
        "results/samples/pe-{sample}/bam/{sample}-clean.sorted.dupsmarked.bam"
    threads: min(cpu_count(), 2)
    shell:
        "samtools view -f 3 -F 3852 -@ {threads} -b  {input} > {output}"
