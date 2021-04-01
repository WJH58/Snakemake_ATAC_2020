#########################
#Sets of rules for HMMRATAC
#Will only work with Paired-end ATAC-data
#--no-mixed argument is required in bowtie2 rules otherwise HMMRATAC will throw an error message.
#HMMRATAC will not work with the test dataset. The sample size is too small.
#########################



rule download_hmmratac:
    output:
        touch(WORKING_DIR + "hmmratac.java")
    params:
        hmmratac = config["HMMRATAC"]
    log:
        RESULT_DIR + "logs/hmmratac/download_hmmratac_github.log"
    shell:
        "wget {params.hmmratac}"

rule hmmratac:
    input:
        map_sorted          = WORKING_DIR + "sort/{samples}.sorted.bam",
        map_sorted_indexed  = WORKING_DIR + "sort/{samples}.sorted.bam.bai"
    output:
        gapped_peak         = RESULT_DIR + "hmmratac/{samples}_peaks.gappedPeak"
    params:
        genome_info         = str(config['hmmratac']['genome.info']),
        output_preflix      = "{samples}"
    log:
        RESULT_DIR + "logs/hmmratac/{samples}_hmmratac.log"
    shell:
        """
        java -Xmx46G -jar HMMRATAC_V1.2.10_exe.jar -b {input.map_sorted} \
        -i {input.map_sorted_indexed} \
        -g {params.genome_info} \
        -o {params.output_preflix} \
        2>{log}
        """
