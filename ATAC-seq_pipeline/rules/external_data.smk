#########################
#Sets of rules for getting external data
#########################

rule download_genome:
    output:
        WORKING_DIR + "reference",
    message:
        "Downloading reference genome."
    shell:
        "wget -O {output} {GENOME_FASTA_URL}"


rule download_gene_gtf:
    output:
        gft         =       WORKING_DIR + "gtf_gene.gtf"
    message:
        "Downloading annotation file."
    shell:
        "curl {GENE_GTF_URL} | gunzip -c > {output}"

rule index_genome:
    input:
        WORKING_DIR + "reference",
    output:
        [WORKING_DIR + "genome." + str(i) + ".bt2" for i in range(1,4)],
        WORKING_DIR + "genome.rev.1.bt2",
        WORKING_DIR + "genome.rev.2.bt2"
    message:
        "Indexing Reference genome."
    params:
        WORKING_DIR + "genome"
    conda:
        "../envs/samtools_bowtie.yaml"
    shell:
        "bowtie2-build --threads 10 {input} {params}"
