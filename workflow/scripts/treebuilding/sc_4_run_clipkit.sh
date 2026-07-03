#!/bin/bash

# --- SLURM Settings ---
#SBATCH --job-name=clipkit    # Job name
#SBATCH --partition=main          # Partition/Queue name
#SBATCH --nodes=1                     # Run on a single node
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=4             # Number of CPU cores per file
#SBATCH --mem=100G                      # Memory limit
#SBATCH --time=00:30:00               # Time limit (hrs:min:sec)
#SBATCH --output=slurmlog/clipkit.log           # Standard output and error log (%j = JobID)

FN_IN=mafft/fbp1_btuF_desE_hmuT_fecB_btuF2_mtsA_clustered_source_tax_subdb1pct-1esz_2chu_3eiw_2r79_1n2z/fbp1_btuF_desE_hmuT_fecB_btuF2_mtsA_clustered_source_tax_subdb1pct-1esz_2chu_3eiw_2r79_1n2z.faa.aln.known_trim
FN_OUT=${FN_IN}.clipkit

clipkit \
    $FN_IN \
    -o $FN_OUT \
    -m kpi-gappy \
    --gaps 0.9 \
    -l \
    -t $SLURM_CPUS_PER_TASK
