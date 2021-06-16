#########################
#MultiQC rule 
#########################

rule multiQC:
    output:
        multiqc = RESULT_DIR + "multiqc/multiqc_report.html"
    conda:
        "../envs/multiqc.yaml"
    log:
        RESULT_DIR + "logs/multiqc/multiqc.log"
    params:
        out_folder = RESULT_DIR + "multiqc/"
    message:
        "Giving summary of fastqc report across samples."
    shell:
        "multiqc . -o {params.out_folder} 2>{log}"
