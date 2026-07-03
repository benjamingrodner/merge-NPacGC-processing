#!/usr/bin/env bash

# --- SLURM Settings ---
#SBATCH --job-name=merge    # Job name
#SBATCH --partition=main          # Partition/Queue name
#SBATCH --nodes=1                     # Run on a single node
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=10             # Number of CPU cores per file
#SBATCH --mem=10G                      # Memory limit
#SBATCH --time=02:00:00               # Time limit (hrs:min:sec)
#SBATCH --output=slurmlog/merge.log           # Standard output and error log (%j = JobID)

######################
# Merge aaseqs
WRK=/scratch/bgrodner/iron_ko_contigs/sidero_receptors/fbp1/alignment/groups_source_taxa

SUBDB=1
PCT_SUBDB=${SUBDB#*.}
BN=fbp1_desE_clustered_source_tax-1esz_2chu_3eiw_2r79_1n2z
DIR_OUT=mafft/"$BN"
mkdir -p "$DIR_OUT"
FN_MERGE="$DIR_OUT"/"$BN".faa

KNOWNS=( 1esz 2chu 3eiw 2r79 1n2z )
DIRS_GENE=( fbp1 desE )
DIRS_TYP=( Database  Environmental_isolate  Environmental_metatranscriptome )
ident="0.9"
cov="0.5"
id="${ident#*.}"0
c="${cov#*.}"0
DIR_CLUST=mmseqsi"${id}"c"${c}"covmode1

TTNAMES=( Bacillariophyta Haptophyta Chlorophyta Dinophyceae Pelagomonas )
TTPATTERN=$(IFS="|"; echo "${TTNAMES[*]}")

# conda 
source activate clipkit

# Prep file
> "$FN_MERGE"

# Add known sequences
DIR_KNOWN=/scratch/bgrodner/iron_ko_contigs/sidero_receptors/fbp1/structure/known_substrate_binding_proteins/berntesson_et_al_2010/fastas
for known in "${KNOWNS[@]}"; do
    cat ${DIR_KNOWN}/${known}.faa >> $FN_MERGE
done

# Add clustered sequences from env
for dir_gene in "${DIRS_GENE[@]}"; do
    echo -e "\n$dir_gene"
    # if [[ "$dir_gene" == 'hmuT' ]]; then
    #     dirs_typ=( Environmental_isolate  Environmental_metatranscriptome )
    # else
    #     dirs_typ=$DIRS_TYP
    # fi
    for dir_typ in "${DIRS_TYP[@]}"; do
        echo -e "\t$dir_typ"
        for fn_repseq in "$WRK"/"$dir_gene"/"$dir_typ"/"$DIR_CLUST"/*_rep_seq.fasta; do

            # Subset the database sequences
            if [[ "$dir_typ" == 'Database' ]]; then
                # But only if they are not in the target sequences of interest
                if [[ "$fn_repseq" =~ $TTPATTERN ]]; then
                    cat "$fn_repseq" >> "$FN_MERGE"
                    echo -e "\t\t Standard cat ${fn_repseq}"
                else
                    echo "seqkit sample -p $SUBDB $fn_repseq >> $FN_MERGE"
                    seqkit sample -p "$SUBDB" "$fn_repseq" >> "$FN_MERGE"
                    echo -e "\t\t Subset ${fn_repseq}"
                fi
            else
                # echo -e "\t\t ${fn_repseq}"
                cat "$fn_repseq" >> "$FN_MERGE"
                echo -e "\t\t Standard cat ${fn_repseq}"
            fi
        done
    done
done

# # Remove duplicates
# FN_MERGE_RMDUP="$FN_MERGE".rmdup
# seqkit rmdup "$FN_MERGE" -o "$FN_MERGE_RMDUP"

# # Run mafft
# mafft --auto \
#     --thread $SLURM_CPUS_PER_TASK \
#     "$FN_MERGE_RMDUP" \
#     > "$FN_MERGE_RMDUP".aln
