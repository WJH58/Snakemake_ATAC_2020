#########################
#Sets of rules for trimming, mapping, sorting and indexing of ATAC-seq reads.
#########################


rule trimmomatic:
    input:
        reads = get_fastq,
        adapters = config['trimmomatic']["adapters"]
    output:
        forward_reads   = WORKING_DIR + "trimmed/{samples}_forward.fastq.gz",
        reverse_reads   = WORKING_DIR + "trimmed/{samples}_reverse.fastq.gz",
        forwardUnpaired = temp(WORKING_DIR + "trimmed/{samples}_forward_unpaired.fastq.gz"),
        reverseUnpaired = temp(WORKING_DIR + "trimmed/{samples}_reverse_unpaired.fastq.gz")
    log:
        RESULT_DIR + "logs/trimmomatic/{samples}.log"
    params:
        seedMisMatches =            str(config['trimmomatic']['seedMisMatches']),
        palindromeClipTreshold =    str(config['trimmomatic']['palindromeClipTreshold']),
        simpleClipThreshhold =      str(config['trimmomatic']['simpleClipThreshold']),
        LeadMinTrimQual =           str(config['trimmomatic']['LeadMinTrimQual']),
        TrailMinTrimQual =          str(config['trimmomatic']['TrailMinTrimQual']),
        windowSize =                str(config['trimmomatic']['windowSize']),
        avgMinQual =                str(config['trimmomatic']['avgMinQual']),
        minReadLen =                str(config['trimmomatic']['minReadLength']),
        phred = 		            str(config["trimmomatic"]["phred"])
    threads: 10
    conda:
        "../envs/trimmomatic.yaml"
    message:
        "Trimming reads for {wildcards.samples}"
    shell:
        """
        trimmomatic PE \
        -threads 10 \
        -phred33 \
         {input.reads} \
         {output.forward_reads} {output.forwardUnpaired} {output.reverse_reads} {output.reverseUnpaired} \
         ILLUMINACLIP:{input.adapters}:{params.seedMisMatches}:{params.palindromeClipTreshold}:{params.simpleClipThreshhold} \
         LEADING:3 \
         TRAILING:3 \
         SLIDINGWINDOW:4:15 \
         MINLEN:40 &>{log}
        """

rule trimmed_fastqc:
    input:
        forward_reads = WORKING_DIR + "trimmed/{samples}_forward.fastq.gz",
        reverse_reads = WORKING_DIR + "trimmed/{samples}_reverse.fastq.gz"
    output:
        RESULT_DIR + "trimmed_fastqc/{samples}_{direction}_fastqc.html"
    params:
        RESULT_DIR + "trimmed_fastqc/"
    conda:
        "../envs/fastqc.yaml"
    message:
        "Quality check after trimming for {wildcards.samples}"
    log:
        RESULT_DIR + "logs/trimmed_fastqc/{samples}_{direction}.fastqc.log"
    shell:
        "fastqc --outdir={params} {input.forward_reads} {input.reverse_reads} &>{log}"

rule fastqc:
    input:
        fwd = expand(DATA_DIR + "sample_{numbers}_R1.fastq.gz", numbers = ['8','12','4_3','4_1']),
        rev = expand(DATA_DIR + "sample_{numbers}_R2.fastq.gz", numbers = ['8','12','4_3','4_1'])
    output:
        expand(RESULT_DIR + "fastqc/sample_{numbers}_{R}_fastqc.html", numbers = ['8','12','4_3','4_1'], R=['R1', 'R2'])
    log:
        expand(RESULT_DIR + "logs/fastqc/{samples}.fastqc.log", samples = SAMPLES)
    params:
        RESULT_DIR + "fastqc/"
    conda:
        "../envs/fastqc.yaml"
    message:
        "Quality check for {wildcards.samples}"
    shell:
        "fastqc --outdir={params} {input.fwd} {input.rev} &>{log}"

rule mapping:
    input:
        forward_reads = WORKING_DIR + "trimmed/{samples}_forward.fastq.gz",
        reverse_reads = WORKING_DIR + "trimmed/{samples}_reverse.fastq.gz",
        forwardUnpaired = WORKING_DIR + "trimmed/{samples}_forward_unpaired.fastq.gz",
        reverseUnpaired = WORKING_DIR + "trimmed/{samples}_reverse_unpaired.fastq.gz",
        index = [WORKING_DIR + "genome." + str(i) + ".bt2" for i in range(1,4)]
    output:
        mapped = WORKING_DIR + "mapped/{samples}.bam",
        unmapped = [WORKING_DIR + "unmapped/{samples}.fq." + str(i) +".gz" for i in range(1,2)],
    params:
        bowtie          = " ".join(config["bowtie2"]["params"].values()), #take argument separated as a list separated with a space
        index           = WORKING_DIR + "genome",
        unmapped        = WORKING_DIR + "unmapped/{samples}.fq.gz"
    threads: 10
    conda:
        "../envs/samtools_bowtie.yaml"
    message:
        "Mapping {wildcards.samples} to reference genome."
    log:
        RESULT_DIR + "logs/bowtie/{samples}.log"
    shell:
        """
        bowtie2 {params.bowtie} --threads {threads} \
        -x {params.index} -1 {input.forward_reads} \
        -2 {input.reverse_reads} \
        -U {input.forwardUnpaired},{input.reverseUnpaired} \
        --un-conc-gz {params.unmapped} --no-mixed | samtools view -Sb - > {output.mapped} 2>{log}
        """

rule sort:
    input:
        mapped  = WORKING_DIR + "mapped/{samples}.bam"
    output:
        sorted  = WORKING_DIR + "sort/{samples}.sorted.bam",
        bai     = WORKING_DIR + "sort/{samples}.sorted.bam.bai"
    conda:
        "../envs/samtools_bowtie.yaml"
    message:
        "Sorting {wildcards.samples} bam file."
    log:
        RESULT_DIR + "logs/sort/{samples}.log"
    shell:
        """
        samtools sort {input} -o {output.sorted} &>{log}
        samtools index {output.sorted}
        """

rule deduplication:
    input:
        map_sorted = WORKING_DIR + "sort/{samples}.sorted.bam"
    output:
        dedup = WORKING_DIR + "dedup/{samples}.dedup.bam",
        stats = WORKING_DIR + "dedup/{samples}.dedup.stats"
    conda:
        "../envs/picard.yaml"
    message:
        "Removing duplicates for {wildcards.samples} sorted bam file."
    log:
        RESULT_DIR + "logs/dedup/{samples}.log"
    shell:
        """
        picard MarkDuplicates -I {input} -O {output.dedup} -M {output.stats} 2>{log}
        """
