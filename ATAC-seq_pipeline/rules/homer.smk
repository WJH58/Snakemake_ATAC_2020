#########################
#Use Homer to annotate files.
#########################
rule annotation:
    input:
        RESULT_DIR + "macs2/{samples}_peaks.narrowPeak"
    output:
        annotation               = RESULT_DIR + "annotation/{samples}_annotate_homer.txt"
    message:
        "Using Homer to annotate {wildcards.samples} peak-calling file"
    log:
        RESULT_DIR + "logs/annotation/{samples}_homer.log"
    conda:
        "../envs/homer.yaml"
    shell:
        "annotatePeaks.pl {input} mm10 > {output}"
