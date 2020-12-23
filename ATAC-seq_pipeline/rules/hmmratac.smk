rule get_genome_info:
    input:
        map_sorted = WORKING_DIR + "sort/{samples}.sorted.bam"
    output:
        genome_info = WORKING_DIR + "genome_info/{samples}.genome.info"
    conda:
        "../envs/perl.yaml"
    shell:
        """
        samtools view -H {input} | grep SQ | cut -f 2-3 | cut -d ':' -f 2 | cut -f 1 > tmp
        samtools view -H {input} | grep SQ | cut -f 2-3 | cut -d ':' -f 2,3 | cut -d ':' -f 2 > tmp2
        paste tmp tmp2 > {output}
        rm tmp tmp2
        """

rule peak_calling:
    input:
        map_sorted = WORKING_DIR + "sort/{samples}.sorted.bam",
        map_sorted_indexed = WORKING_DIR + "index/{samples}.sorted.bam.bai"
    output:
        gapped_peak = RESULT_DIR + "peaks/{samples}_peaks.gappedPeak",
        Name_log = RESULT_DIR + "peaks/{samples}.log"
    params:
        genome_info = str(config['hmmratac']['genome.info']),
        output_preflix = "{samples}"
    log:
        RESULT_DIR + "logs/peak_calling/{samples}.log"
    conda:
        "../envs/hmmratac.yaml"
    shell:
        """
        HMMRATAC -b {input.map_sorted} -i {input.map_sorted_indexed} -g {params.genome_info} -o {params.output_preflix}
        """
