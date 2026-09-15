#!/usr/bin/env bash

# Set the job name
#SBATCH --job-name=snakemake_master

# Request one core for the master process
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

# Set the maximum time the master job can run (e.g., 1 hour, or until the workflow finishes)
#SBATCH --time=12:00:00

# Set the partition/queue to use
#SBATCH --partition=main

# Set the memory for the master process
#SBATCH --mem=300G

# Redirect the master job's stdout and stderr to a file
#SBATCH --output=slurm_logs/snakemake_master.%j.out
#SBATCH --error=slurm_logs/snakemake_master.%j.err

IMAGE=docker://benjamingrodner/get_metat_dicts
SNAKEFILE=workflow/Snakefile
apptainer run \
    --no-home \
    --bind /mnt/nfs/projects/armbrust-metat \
    --bind /scratch/bgrodner \
    $IMAGE \
    snakemake \
        --snakefile $SNAKEFILE \
        --configfile config.yaml \
        --jobs 32 \
        --rerun-triggers mtime \
        --rerun-incomplete \
        -p
        # -R results/tmp/tax_lineage/g1-st-am-pa-tax_lineage.parquet \
        # -R results/tmp/merge_counts/g2-st-am-pa-merge_counts.parquet \
        # --jobs $SLURM_CPUS_PER_TASK \
        # --unlock \

        # --latency-wait 60 \
        # --resources mem_gb=200 \
        # --jobs 1000 \
        # --cluster-config cluster_config.yaml \
        # --cluster "sbatch --partition={cluster.partition} --time={cluster.time} --mem={cluster.mem} --cpus-per-task={cluster.n} --output=slurm_logs/{rule}-%j.out --error=slurm_logs/{rule}-%j.err" \

