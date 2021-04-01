#! /bin/bash

#SBATCH --partition=highmem
#SBATCH  --job-name=mini_snakemake
#SBATCH --mail-type=ALL
#SBATCH --mail-user j.wang@lumc.nl
#SBATCH -t 48:00:00
#SBATCH --mem=48G


echo Start time : `date`
snakemake -p \
        --snakefile Snakefile \
        --latency-wait 60 \
        --wait-for-files \
        --rerun-incomplete \
        --use-conda \
        --cluster "sbatch --parsable --partition=highmem --mem=60g --ntasks=1 --cpus-per-task=8 --time=60:00:00 --hint=multithread" \
        --cluster-status "./slurm-cluster-status.py" \
       	--jobs 30


echo End time : `date`
