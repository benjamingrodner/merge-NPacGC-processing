#!/usr/bin/env bash

# --- SLURM Settings ---
#SBATCH --job-name=build_tree    # Job name
#SBATCH --partition=main          # Partition/Queue name
#SBATCH --nodes=1                     # Run on a single node
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per file
#SBATCH --mem=8G                      # Memory limit
#SBATCH --time=00:30:00               # Time limit (hrs:min:sec)
#SBATCH --output=slurmlog/build_tree_%j.log           # Standard output and error log (%j = JobID)


######################
# Basename
BN=


######################
# Merge aaseqs

FN_OUT=fbp1-known_struct-Chlorophyta-contigs_of_interest-renamed.faa

# Add fbp1
FN_FBP1=/scratch/bgrodner/iron_ko_contigs/sidero_receptors/fbp1/FBP1_analogs/Phatr3_J46929.p1.fasta
cat $FN_FBP1 > $FN_OUT

# Add known sequences
DIR_KNOWN=/scratch/bgrodner/iron_ko_contigs/sidero_receptors/fbp1/structure/known_substrate_binding_proteins/berntesson_et_al_2010/fastas
KNOWNS=('1esz' '2chu' '3eiw')
for known in "${KNOWNS[@]}"; do
    cat ${DIR_KNOWN}/${known}.faa >> $FN_OUT
done

# add contigs of interest
FN_IN=/scratch/bgrodner/iron_ko_contigs/sidero_receptors/fbp1/alignment/aaseqs_merged/Chlorophyta-contigs_of_interest-renamed.faa
cat $FN_IN >> $FN_OUT    

echo -e "Wrote aaseqs\n\n\n"
# DIR_CRYSTAL=/scratch/bgrodner/iron_ko_contigs/sidero_receptors/fbp1/structure/known_substrate_binding_proteins/berntesson_et_al_2010/fastas/rename
# for fn_fasta in ${DIR_CRYSTAL}/*.faa; do
#     cat $fn_fasta >> $FN_OUT
# done





######################
# run mafft

echo -e "Running mafft...\n\n\n"

FN_IN=$FN_OUT
FN_OUT=${FN_OUT%.faa}.aln
mafft --auto \
    $FN_IN \
    > $FN_OUT

echo -e "Ran mafft\n\n\n"


######################
# run clipkit

echo -e "Running clipkit...\n\n\n"


FN_IN=$FN_OUT
FN_OUT=${FN_IN}.clipkit

clipkit \
    $FN_IN \
    -o $FN_OUT \
    -l \
    -m kpi-gappy \
    --gaps 0.75 

echo -e "Ran clipkit\n\n\n"


######################
# filter_alignment

MINLEN=80
python sc_filter_alignment.py $FN_OUT $MINLEN

echo -e "Filtered alignments\n\n\n"


######################
# run fasttree

echo -e "Running fasttree...\n\n\n"


FN_IN=$FN_OUT.minlen${MINLEN}
FN_OUT=${FN_IN}.fasttree

fasttreeMP \
    -fastest \
    $FN_IN \
    > $FN_OUT


echo -e "Ran fasttree\n\n\n"
