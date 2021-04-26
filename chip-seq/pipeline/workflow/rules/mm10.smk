from multiprocessing import cpu_count


rule download_mm10_chrominfo:
    output:
        "resources/mm10/chrominfo.txt"
    shell:
        "wget -O {output} https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes"


rule download_mm10_genome:
    output:
        "resources/mm10/genome.fa"
    shell:
        "workflow/scripts/download-genome.sh resources/mm10"


rule build_mm10_minimap2_index:
    input:
        "resources/mm10/genome.fa"
    output:
        "resources/mm10/genome.mmi"
    threads: cpu_count()
    params:
        mode="sr"
    shell:
        "minimap2 -x {params.mode} -t {threads} -d {output} {input}"


rule download_mm10_blacklist:
    output:
        "resources/mm10/mm10-blacklist.v2.bed.gz"
    shell:
        "wget -O {output} https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm10-blacklist.v2.bed.gz"
