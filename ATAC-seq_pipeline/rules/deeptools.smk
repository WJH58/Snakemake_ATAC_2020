#########################
#Sets of rules for deeptools postprocessing. All parameters can be modified in the configuration file.
#########################

rule bamCoverage:
    input:
        map_sorted              = WORKING_DIR + "sort/{samples}.sorted.bam"
    output:
        coverage_track          = RESULT_DIR + "bamCoverage/{samples}.bw"
    params:
        effectiveGenomeSize     = str(config['bamCoverage']['effectiveGenomeSize']),
        normalization           = str(config['bamCoverage']['normalization']),
        binSize                 = str(config['bamCoverage']['binSize'])
    log:
        RESULT_DIR + "logs/bamCoverage/{samples}.bw.log"
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        bamCoverage -b {input} \
        --effectiveGenomeSize {params.effectiveGenomeSize} \
        --normalizeUsing {params.normalization} \
        --ignoreDuplicates \
        -o {output} \
        --binSize {params.binSize}
        """
#bdg: create bedgraph output files
#format=BAMPE: only use properly-paired read alignments


# rule multibigwigSummary need to be expanded to use more than 1 bigwig file
# Right now the file only use 1 file, which is not useful for a comparison
rule multibigwigSummary:
    input:
        lambda wildcards: expand(RESULT_DIR + "bamCoverage/{sample}.bw", sample = SAMPLES)
    output:
        bigwigsummary           = RESULT_DIR + "bigwigsummary/multiBigwigSummary.npz"
    log:
        RESULT_DIR + "logs/bigwigsummary/multiBigwigSummary.log"
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        multiBigwigSummary bins -b {input} \
        -o {output}
        """
#2>{log}
# PlotPCA does not work with a bigwigsummary made of only one sample
rule PCA:
    input:
        bigwigsummary           = RESULT_DIR + "bigwigsummary/multiBigwigSummary.npz"
    output:
        PCAplot                 = RESULT_DIR + "PCA/PCA_PLOT.pdf"
    log:
        RESULT_DIR + "logs/bigwigsummary/PCA_PLOT.log"
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        plotPCA -in {input} \
        -o {output} \
        2>{log}
        """
rule ComputeMatrix:
    input:
        lambda wildcards: expand(RESULT_DIR + "bamCoverage/{sample}.bw", sample = SAMPLES)
    output:
        RESULT_DIR + "computematrix/ComputeMatrix.gz"
    log:
        RESULT_DIR + "logs/computematrix/matrix.log"
    conda:
        "envs/deeptools.yaml"
    params:
        GTF = WORKING_DIR + "gtf_gene.gtf"
    shell:
        "computeMatrix reference-point -S {input} -R {params.GTF} -o {output} -a 3000 -b 3000 2>{log}"

rule plotHeatmap:
    input:
        RESULT_DIR + "computematrix/ComputeMatrix.gz"
    output:
        RESULT_DIR + "heatmap/heatmap_reference_point_genes.pdf"
    conda:
        "envs/deeptools.yaml"
    shell:
        "plotHeatmap -m {input} -o {output}"
