#########################
#Sets of rules for deeptools postprocessing. All parameters can be modified in the configuration file.
#########################
rule bamPEFragmentSize:
    input:
        lambda wildcards: expand(WORKING_DIR + "mapped/{sample}.bam", sample = SAMPLES)
    output:
        histogram               = RESULT_DIR + "bamPEFragmentSize/histogram.png"
    conda:
        "../envs/deeptools.yaml"
    message:
        "Calculating fragment sizes for read pairs."
    log:
        RESULT_DIR + "logs/bamPEFragmentSize/bamPEFragmentSize.log"
    params:
        numberOfProcessors      = str(config['bamPEFragmentSize']['numberOfProcessors'])
        binSize                 = str(config['bamPEFragmentSize']['binSize'])
        plotFileFormat          = str(config['bamPEFragmentSize']['plotFileFormat'])
        samplesLabel            = str(config['bamPEFragmentSize']['samplesLabel'])
        plotTitle               = str(config['bamPEFragmentSize']['plotTitle'])
    shell:
        """
        bamPEFragmentSize -b {input} \
        -o {output} \
        --numberOfProcessors {params.numberOfProcessors} \
        --binSize 1000 {params.binSize} \
        --plotFileFormat {params.plotFileFormat} \
        --samplesLabel {params.samplesLabel} \
        -T {params.plotTitle} \
        --table 2>{log}
        """

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
    message:
        "Converting {wildcards.samples} into bigwig files."
    shell:
        """
        bamCoverage -b {input} \
        --effectiveGenomeSize {params.effectiveGenomeSize} \
        --normalizeUsing {params.normalization} \
        --ignoreDuplicates \
        -o {output} \
        --binSize {params.binSize} 2>{log}
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
    message:
        "Computing the read coverage into a numpy array."
    shell:
        """
        multiBigwigSummary bins -b {input} \
        -o {output}
        """

rule plotCorrelation:
    input:
        bigwigsummary           = RESULT_DIR + "bigwigsummary/multiBigwigSummary.npz"
    output:
        correlation             = RESULT_DIR + "correlation/correlation.pdf"
    log:
        RESULT_DIR + "logs/correlation/correlation.log"
    params:
        corMethod               = str(config['plotCorrelation']['corMethod'])
        whatToPlot              = str(config['plotCorrelation']['whatToPlot'])
    message:
        "Plotting correlation."
    shell:
        """
        plotCorrelation --corData {input} \
        --corMethod {params.corMethod} \
        --whatToPlot {params.whatToPlot} \
        -o {output} \
        --skipZeros 2>{log}
        """

# 2>{log}
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
    message:
        "Plotting PCA"
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
        matrix                   = RESULT_DIR + "computematrix/ComputeMatrix.gz"
    log:
        RESULT_DIR + "logs/computematrix/matrix.log"
    conda:
        "../envs/deeptools.yaml"
    message:
        "Computing matrix for samples with {params.binSize} windows and {params.afterRegionStartLength}bp around {params.referencePoint}."
    params:
        GTF                      = WORKING_DIR + "gtf_gene.gtf"
        binSize                  = str(config['ComputeMatrix']['binSize'])
        afterRegionStartLength   = str(config['ComputeMatrix']['afterRegionStartLength'])
        beforeRegionStartLength  = str(config['ComputeMatrix']['beforeRegionStartLength'])
        referencePoint           = str(config['ComputeMatrix']['referencePoint'])
        numberOfProcessors       = str(config['ComputeMatrix']['numberOfProcessors'])
    shell:
        """
        computeMatrix reference-point \
        -S {input} -R {params.GTF} -o {output} \
        -a {params.afterRegionStartLength} \
        -b {params.beforeRegionStartLength} \
        --referencePoint {params.referencePoint} \
        --binSize {params.binSize} \
        --skipZeros \
        2>{log}
        """

rule plotHeatmap:
    input:
        RESULT_DIR + "computematrix/ComputeMatrix.gz"
    output:
        heatmap                  = RESULT_DIR + "heatmap/heatmap_reference_point_genes.pdf"
    log:
        RESULT_DIR + "logs/heatmap/heatmap.log"
    conda:
        "../envs/deeptools.yaml"
    message:
        "Plotting heatmap."
    params:
        dpi                      = str(config['plotHeatmap']['dpi'])
        yMin                     = str(config['plotHeatmap']['yMin'])
        yMax                     = str(config['plotHeatmap']['yMax'])
        refPointLabel            = str(config['plotHeatmap']['refPointLabel'])
        colorList                = str(config['plotHeatmap']['colorList'])
    shell:
        """
        plotHeatmap -m {input} -o {output} \
        --yMin {params.yMin} \
        --yMax {params.yMax} \
        --refPointLabel {params.refPointLabel} \
        --colorList {params.colorList} 2>{log}"
        """

rule plotProfile:
    input:
        RESULT_DIR + "computematrix/ComputeMatrix.gz"
    output:
        profile                  = RESULT_DIR + "profile/profile_reference_point_genes.pdf"
    log:
        RESULT_DIR + "logs/profile/profile.log"
    conda:
        "../envs/deeptools.yaml"
    message:
        "Plotting profile."
    params:
        dpi                      = str(config['plotProfile']['dpi'])
        colors                   = str(config['plotProfile']['colors'])
    shell:
        """
        plotProfile -m {input} \
        -out {output} \
        --dpi {params.dpi} \
        --colors {params.colors} \
        --perGroup 2>{log}
        """
