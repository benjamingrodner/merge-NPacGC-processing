#!/usr/bin/env bash

# --- SLURM DIRECTIVES FOR THE MASTER SNAKEMAKE JOB ---
# This job runs the Snakemake orchestrator, not the individual rules.

# Set the job name
#SBATCH --job-name=snakemake_master

# Request one core for the master process
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10

# Set the maximum time the master job can run (e.g., 1 hour, or until the workflow finishes)
#SBATCH --time=01:00:00

# Set the partition/queue to use
#SBATCH --partition=main

# Set the memory for the master process
#SBATCH --mem=400G

# Redirect the master job's stdout and stderr to a file
#SBATCH --output=slurm_logs/snakemake_master.%j.out
#SBATCH --error=slurm_logs/snakemake_master.%j.err

IMAGE=docker://benjamingrodner/get_metat_dicts
SNAKEFILE=/scratch/bgrodner/repo-armbrust-metat-search/Snakefile_big_tables
apptainer run \
    --no-home \
    --bind /mnt/nfs/projects/armbrust-metat \
    --bind /scratch/bgrodner \
    $IMAGE \
    snakemake \
        --snakefile $SNAKEFILE \
        --configfile config.yaml \
        --jobs $SLURM_CPUS_PER_TASK \
        --rerun-incomplete \
        -p
        # --unlock \
        # --rerun-triggers mtime \

        # --latency-wait 60 \
        # --resources mem_gb=200 \
        # --jobs 1000 \
        # --cluster-config cluster_config.yaml \
        # --cluster "sbatch --partition={cluster.partition} --time={cluster.time} --mem={cluster.mem} --cpus-per-task={cluster.n} --output=slurm_logs/{rule}-%j.out --error=slurm_logs/{rule}-%j.err" \

