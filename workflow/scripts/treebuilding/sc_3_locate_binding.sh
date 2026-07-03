#!/bin/bash

source activate clipkit

FN_ALN=mafft/fbp1_btuF_desE_hmuT_fecB_btuF2_mtsA_clustered_source_tax_subdb1pct-1esz_2chu_3eiw_2r79_1n2z/fbp1_btuF_desE_hmuT_fecB_btuF2_mtsA_clustered_source_tax_subdb1pct-1esz_2chu_3eiw_2r79_1n2z.faa.aln

KNOWNS=( 1esz 2chu 3eiw 2r79 1n2z )

RE_ID="[0-9a-zA-Z]*chain_A"

FN_RESID_SEARCH=known_structs_resids_search.txt

FN_LOC_TABLE="$FN_ALN".known_struct_locs.tab
FN_STARTEND="$FN_ALN".known_struct_startend.csv
FN_TRIM="$FN_ALN".known_trim

# Get start and end positions
seqkit grep -r -p "$RE_ID" "$FN_ALN"  \
    | seqkit replace -p "\s.+" -r "" \
    | seqkit fx2tab \
    | awk '{ 
        match($2, /[^-]/); 
        s_pos=RSTART; 
        match($2, /.*[^-]/); 
        e_pos=RLENGTH; 
        if(min=="" || s_pos < min) min=s_pos; 
        if(e_pos > max) max=e_pos
    } END{print min "," max}' \
    > "$FN_STARTEND"

# Trim alignment file

    
cat "$FN_STARTEND" \
    | (IFS=',' read -r start end; 
        seqkit subseq -r "$start":"$end" "$FN_ALN" > "$FN_TRIM"
    )

# extra
# Get table of locations of left and right boundaries on binding pocket for different known structures
seqkit grep -r -p "$RE_ID" "$FN_ALN" \
    | seqkit locate -t unlimit -r -f "$FN_RESID_SEARCH" \
    > "$FN_LOC_TABLE"

# # leftmost and rightmost 
# awk -F'\t' 'NR > 1 {
#     if(min=="" || $5 < min) min=$5; 
#     if($6 > max) max=$6
# } END{print min,max}' "$FN_LOC_TABLE"

