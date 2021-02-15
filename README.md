# Snakemake pipepine for (Illumina) paired-end ATAC-Seq data

## Overview

This repository contains a paired-end ATAC-seq analysis pipeline. Starting from raw data in fasta (or fastq) format, the pipeline will perform quality control, alignment peak calling and various downstream analysis (including heatmaps, PCA, etc)

The expected outputs of this pipeline are:
**Here make a list with bullet point of the outputs**
The results are a set of files for quality control report, alignment, peak information, consensus and differential peak information.
**Add information about the type of computer required for the pipeline to run (we have a script to launch it on SLURM but it can also work on any computer)**

![DAG](QC_trimmed/dag.png)

## Quick Start


***The actual fist step is to git clone this repository, please add the information about this**

### 1. Conda activate a environment for running Snakemake

This pipeline relies on [conda](https://docs.conda.io/en/latest/miniconda.html) to install packages. First you need to create a conda environment called "snakemake" with the following command: ```conda env create snakemake```, then activate this environment: ```conda activate snakemake```. Run the pipeline in this environment. **You need first to install snakemake via conda in this environment**

The snakefile contains information about the required packages for the pipeline and will ensure the reproducibility of the analysis.

## Content of the repository   **I move this hear**
* The **Snakefile** is the core of the Snakemake workflow. It ensures the reproducibility of the analysis byt setting a set of rules and environments that will produce the desired outputs. This file is hard coded  and should **not** be modified by unexperienced users.

* **Config.yaml** contains all the parameters you want to refer to, such as directory paths, genome file link, rule parameters, etc. You can also save those rarely-modified files such as ```genome.info``` in the configuration file. **Note: all paths in snakemake files are relative!** You can modify Config.yaml to make your pipeline more flexible: genome fasta, gene annotation, deeptools parameters, etc.

* **Envs/** contains conda environment sheets. Every rule that requires tool installation from conda has a separate *.yaml file. ```channels``` refers to which channel of conda you use, and ```dependecies``` refers to the exact version of this tool. To avoid tool conflicts with python, we usually use an older version of a tool. To ensure reproducibility, the environments should not be modified.

* **units.tsv** is a _tab separated value files_ containing information about the experiment name and path to the fastq files relative to the Snakefile. Change this file according to your samples. Fill in the sample names and relative paths when you analyze your own data. **Add information about the importance of the format of the table, it always need to contains the columns ..., ..., ..., and they always need to be tab separated value.** **Add the errors message you get if this file is not in the good format** **File format= fasta/fastq or fasta.gz/fastq.gz.**


### 2. Test it on your computer


**this part misses a lot of informations. You need to add that the repository includes a set of files that can be used to test if the pipeline will run correctly on the computer. Those files are small fastq files of a mouse genomic data. Indicate also how long the test would be on slurm. Indicate also how to run the test and what is the expected outcome**

* **Change genome file** Change 'genome_fasta_url' item in configfile to download the genome that you want to reference.

* **Set binSize value** We set the default binSize value at 1000bp for testing. Remember to change it to small values (e.g.10) when using your data by altering "--binSize" argument of rule _computeMatrix_.

* **Select normalization method** Normalization methods include RPKM, CPM, BPM and RPGC. Our pipeline uses RPKM by default. You can set it after "--normalizeUsing" argument of rule _bamCoverage_.

* **Set matrix parameters** By default the pipeline uses reference-point mode. Alternatively you can change it to scale-regions. Correspondingly you should specify "--regionBodyLength" argument of rule _computeMatrix_.

* **Change units.tsv** Change the sample names and relative paths accordingly. The pipeline includes samples, so you can use it to quickly test whether all environments are installed.

**How to make a test run**

**Output of a test run**

### 3. Snakemake execution

Snakemake workflow management system is a tool to create reproducible data analysis. It executes steps by reading _Snakefile_. Install Snakemake [here](https://snakemake.readthedocs.io/en/stable/).
* Test run
After activating snakemake environment, use command ```snakemake -np``` to perform a dry run. If there is no error messages, you can start a real run.
* Real run
Within the folder containing the Snakefile, simply run this command line ```Snakemake --use-conda --cores n ```. **n is the number of cores you want to dedicate for the analysis**

**Command lines to run also on slurm**

## Output description

 The desired output of this pipeline are:

* **Quality check files of raw data and trimmed data**. Reads are processed and quality reports are generated by fastqc tool. Fastqc reports are in html format. **indicate folder**

* **Alignment files along with indexes**. BAM is short for "Binary Alignment Map", which is a compressed binary version of SAM (Sequence Alignment Map) file. BAM files are quite large so one could decide to have them erased after all analysis is done。They can be useful for other downstream analysis (e.g deeptools). **indicate folder**

* **Peaks**. These files gather the information processed by MACS2. _.narrowPeak_ files includes standard BED files and statistical significance information. **indicate folder**

* **Coverage tracks**. _.bigwig_ files are used for visualizetion on Genome Browser. Including: _NARROWPEAK_, _COVERAGE_TRACK_, _BIGWIGSUMMARY_. **indicate folder**

* **Deeptools visualization of reads over genomic features** (e.g. Promoters, TSS, intergenic regions, etc.). Including: _PCAPLOT_, _HEATMAP_. **indicate folder**


## What do rules do?
* pre-processing.smk:
  * **rule trimmomatic:** Trim the adapters and N bases.
  * **rule trimmed_fastqc:** As confirmation, check the quality again after reads are trimmed.
  * **rule fastqc:** Check the quality of sequenced reads.
  * **rule mapping:** Align trimmed reads with indexed genome.
  * **rule sort&index:** Sort the alignments of BAM file based on chromosomes and index them.
  * **rule deduplication:** Mark the duplicates and remove them.

* macs2.smk:
  * **rule call_narrow_peaks:** Identify reads-enriched regions, namely accessible regions using Macs2. The output files include _.narrowPeak_, which is used for downstream analysis (e.g. deeptools).

* external_data.smk:
  * **rule download_genome:** Download reference genome.
  * **rule download_gene_gtf:** Download gene annotations file.
  * **rule index_genome:** Index reference genome to allow the aligner to narrow down the potential origin of a query sequence within the genome. This saves both time and memory.

* deeptools.smk:
  * **rule bamCoverage:** Convert BAM files to bigwig files that can be visualized on genome browser.
  * **rule multibigwigSummary:** Compute the average scores for each of the files in every genomic region. The output could be used by plotPCA for visualization.
  * **rule PCA:** The PCA plot shows whether samples display greater variability between experimental conditions than between replicates of the same treatment.
  * **rule ComputeMatrix:** Calculate scores per genome regions and prepares an intermediate file that can be used by plotHeatmap or plotProfile.
  * **rule plotHeatmap:** Visualization of peaks over reference point (by default: TSS).

## Configuration file
* working_dir: a directory that contains temporary files.
* result_dir: a directory taht contains desired output files
* data_dir: a directory that stores all the sample data. Sample paths in units.tsv direct here
* units.tsv: a table that indicates sample names and paths. Paired-end sequencing data are accepted and forward reads should be in the column "fq1" and reverse reads should be in the column "fq2". This table should be **tab delimited**
* genome_fasta_url: the link to reference genome
* adapters: adapters.fasta file that is used by trimmomatic. This file is fixed as long as the reads are sequenced by Illumina.
