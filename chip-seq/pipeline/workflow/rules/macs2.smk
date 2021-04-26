####################################################################################
# predictd
####################################################################################

def _macs2_predictd_params(wildcards):
    return {"gsize": samples.organism(wildcards.sample)}


rule macs2_predictd:
    input:
        "results/samples/se-{sample}/bam/{sample}-clean.sorted.dupsmarked.bam"
    output:
        "results/samples/se-{sample}/qc/macs2-fragment-size.txt"
    params:
        p = _macs2_predictd_params
    shell:
        "macs2 predictd --gsize={params.p[gsize]} --rfile=/dev/null -i {input} &> {output}"

####################################################################################
# callpeak
####################################################################################

def _macs2_callpeak_input(wildcards):
    return Design(wildcards.design).dependencies('clean')



rule macs2_callpeak:
    input: unpack(_macs2_callpeak_input)
    output:
        dir=directory("results/peaks/{design}-macs2"),
        peaks="results/peaks/{design}-macs2.narrowPeak"
    threads: 1
    run:
        flags = peak_calling.macs2.flags(wildcards.design)
        
        shell(
            "macs2 callpeak {flags} --outdir={output.dir} --name {wildcards.design}-macs2 "
            "-t {input.treatment} -c {input.control}"
        )
        shell("cp {output.dir}/{wildcards.design}-macs2_peaks.narrowPeak {output.peaks}")

####################################################################################
# signal
####################################################################################

def _macs2_signal_input(wildcards):
    deps = _macs2_callpeak_input(wildcards)
    deps["chrominfo"] = "resources/mm10/chrominfo.txt"
    return deps


rule macs2_signal:
    input: unpack(_macs2_signal_input)
    output:
        fe="results/signal/{design}-macs2-fe.bw",
        control="results/signal/{design}-macs2-control.bw",
        treatment="results/signal/{design}-macs2-treatment.bw",
        tmpdir=temp(directory("/tmp/macs2-signal-{design}/"))
    params:
        name="tmp"
    threads: 1
    run:
        flags = peak_calling.macs2.flags(wildcards.design)
        shell(
            "macs2 callpeak -t {input.treatment} -c {input.control} " 
            "{flags} --bdg --SPMR --name={params.name} --outdir={output.tmpdir}"
        )
        treatment = f"{output.tmpdir}/{params.name}_treat_pileup.bdg"
        control = f"{output.tmpdir}/{params.name}_control_lambda.bdg"
        fe = f"{output.tmpdir}/{params.name}_fe.bdg"
        shell("macs2 bdgcmp -m FE -t {treatment} -c {control} -o {fe}")
        
        for bdg, saveto in [(treatment, output.treatment), (control, output.control), (fe, output.fe)]:
            clipped = f"{output.tmpdir}/{params.name}_tmp_clipped.bdg"
            shell(
                "bedtools slop -i {bdg} -g {input.chrominfo} -b 0 | "
                "bedClip stdin {input.chrominfo} /dev/stdout | "
                "LC_COLLATE=C sort -k1,1 -k2,2n /dev/stdin > {clipped}"
            )
            shell("bedGraphToBigWig {clipped} {input.chrominfo} {saveto}")
