#!/bin/bash

# --- SLURM Settings ---
#SBATCH --job-name=mmseqs    # Job name
#SBATCH --partition=main          # Partition/Queue name
#SBATCH --nodes=1                     # Run on a single node
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=4             # Number of CPU cores per file
#SBATCH --mem=50G                      # Memory limit
#SBATCH --time=01:00:00               # Time limit (hrs:min:sec)
#SBATCH --output=slurmlog/mmseqs.log           # Standard output and error log (%j = JobID)

set -e
module load mmseqs2
WRK=/scratch/bgrodner/iron_ko_contigs/sidero_receptors/fbp1/alignment/groups_source_taxa



ident="0.9"
cov="0.5"

id="${ident#*.}"0
c="${cov#*.}"0

# DIRS_GENE=( fbp1 btuF )
DIRS_GENE=( fbp1 btuF desE hmuT fecB btuF2 mtsA )
# DIRS_TYP=( Environmental_metatranscriptome )
DIRS_TYP=( Database  Environmental_isolate  Environmental_metatranscriptome )
# EXTS_TYP=( merge_cruises.faa )
EXTS_TYP=( .faa .faa merge_cruises.faa )
# count=0
for dir_gene in "${DIRS_GENE[@]}"; do
    echo -e "\n$dir_gene"
    for j in "${!DIRS_TYP[@]}"; do
        echo -e "\t${DIRS_TYP[$j]}"
        for fn_fasta in "$WRK"/"$dir_gene"/"${DIRS_TYP[$j]}"/*"${EXTS_TYP[$j]}"; do
            # [[ $count -eq 3 ]] && break
            # ((count++))
            dir_fasta=$(dirname "$fn_fasta")
            bn_fasta=$(basename "$fn_fasta")
            dir_clust="$dir_fasta"/mmseqsi"${id}"c"${c}"covmode1
            mkdir -p "$dir_clust"
            bn_clust="$dir_clust"/"${bn_fasta%.faa}"
            fn_rs="$bn_clust"_rep_seq.fasta
            if [[ -f "$fn_fasta" ]]; then # Run if there are fastas in the dir
            # if [[ -f "$fn_fasta" && ! -f "$fn_rs" ]]; then # Run if not run before
                mmseqs easy-cluster \
                    "${fn_fasta}" \
                    "${bn_clust}" \
                    tmp \
                    --min-seq-id "${ident}" \
                    -c "${cov}" \
                    --cov-mode 1 \
                    --cluster-mode 2 
            fi
        done
    done
done

