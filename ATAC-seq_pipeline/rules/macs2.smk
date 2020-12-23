#########################
#Sets of rules for MACS2 peak calling
#########################

rule call_narrow_peaks:
    input:
        map_sorted      = WORKING_DIR + "sort/{samples}.sorted.bam"
    output:
        narrowPeak      = RESULT_DIR + "macs2/{samples}_peaks.narrowPeak"
    params:
        name            = "{samples}",
        format          = str(config['macs2']['format']),
        genomesize      = str(config['macs2']['genomesize']),
        outdir          = str(config['macs2']['outdir'])
    log:
        RESULT_DIR + "logs/macs2/{samples}_peaks.narrowPeak.log"
    conda:
        "../envs/macs2.yaml"
    shell:
        """
        macs2 callpeak -t {input} {params.format} {params.genomesize} \
        --name {params.name} --nomodel --bdg -q 0.05 \
        --outdir {params.outdir}/ 2>{log}
        """
