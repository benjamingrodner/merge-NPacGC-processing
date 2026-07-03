#!/usr/bin/env bash

# --- SLURM Settings ---
#SBATCH --job-name=mafft    # Job name
#SBATCH --partition=main          # Partition/Queue name
#SBATCH --nodes=1                     # Run on a single node
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=10             # Number of CPU cores per file
#SBATCH --mem=100G                      # Memory limit
#SBATCH --time=02:00:00               # Time limit (hrs:min:sec)
#SBATCH --output=slurmlog/mafft_%j.log           # Standard output and error log (%j = JobID)

FN_MERGE=mafft/fbp1_desE_clustered_source_tax-1esz_2chu_3eiw_2r79_1n2z/fbp1_desE_clustered_source_tax-1esz_2chu_3eiw_2r79_1n2z.faa

# Run mafft
mafft --auto \
    --thread $SLURM_CPUS_PER_TASK \
    "$FN_MERGE" \
    > "$FN_MERGE".aln