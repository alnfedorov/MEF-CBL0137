from multiprocessing import cpu_count


rule fastqc:
    input:
        "results/samples/{sample}/fq.gz/{fq}.fq.gz"
    output:
        "results/samples/{sample}/qc/fastqc/{fq}_fastqc.html",
        "results/samples/{sample}/qc/fastqc/{fq}_fastqc.zip"
    threads: 1
    shell:
        "fastqc -t {threads} -o results/samples/{wildcards.sample}/qc/fastqc {input}"


rule samtools_stats:
    input:
        "results/samples/{sample}/bam/{bam}.sorted.dupsmarked.bam"
    output:
        "results/samples/{sample}/qc/samtools-stats/{bam}.stats"
    threads: 1
    shell:
        "samtools stats -@ {threads} {input} > {output}"


rule samtools_flagstat:
    input:
        "results/samples/{sample}/bam/{bam}.sorted.dupsmarked.bam"
    output:
        "results/samples/{sample}/qc/samtools-flagstat/{bam}.flagstat"
    threads: cpu_count()
    shell:
        "samtools flagstat -@ {threads} {input} > {output}"


rule fragment_size_estimation:
    input:
        bam="results/samples/pe-{sample}/bam/{sample}-clean.sorted.dupsmarked.bam",
        bai="results/samples/pe-{sample}/bam/{sample}-clean.sorted.dupsmarked.bam.bai"
    output:
        dist="results/samples/pe-{sample}/qc/bamPEFragmentSize/bamPEFragmentSize.RawFragmentLengths",
        hist="results/samples/pe-{sample}/qc/bamPEFragmentSize/histogram.png"
    threads: cpu_count()
    params:
        label="{sample}"
    shell:
        "bamPEFragmentSize -p {threads} --samplesLabel {params.label} --outRawFragmentLengths {output.dist} "
        "--histogram {output.hist} --bamfiles {input.bam}"

rule multiqc_sample_pe:
    input:
        "results/samples/pe-{sample}/fq.gz/{sample}_1.fq.gz_trimming_report.txt",
        "results/samples/pe-{sample}/fq.gz/{sample}_2.fq.gz_trimming_report.txt",
        "results/samples/pe-{sample}/qc/samtools-stats/{sample}-clean.stats",
        "results/samples/pe-{sample}/qc/samtools-flagstat/{sample}-clean.flagstat", 
        "results/samples/pe-{sample}/qc/fastqc/{sample}_1_fastqc.html",
        "results/samples/pe-{sample}/qc/fastqc/{sample}_2_fastqc.html",
        "results/samples/pe-{sample}/qc/bamPEFragmentSize/bamPEFragmentSize.RawFragmentLengths"
    output:
        report="results/samples/pe-{sample}/qc/multiqc_report.html"
    shell:
        "multiqc -o results/samples/pe-{wildcards.sample}/qc results/samples/pe-{wildcards.sample}"
