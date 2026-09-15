rule get_aaseq_parquets:
    input:
        fns_seqs = lambda w: get_fnseqs(w.batch, input_table)
    output:
        fn_seqs_parq = fmt_seqs_parq
    threads: 8
    resources:
        mem="300"
    shell:
        """
        # From chris' pipeline
        # Write AA ORFs contig name, longest AA ORF name, AA seq to Parquet file
        echo "$(date): Writing AA contig name,  AA name, AA seq to parquet file {output.fn_seqs_parq}"
        AA_ORFS_SEQS_SQL=$(mktemp --suffix=.sql)
        cat > "$AA_ORFS_SEQS_SQL" <<EOF
        COPY (
            SELECT
                regexp_replace(split_part(Column0, ' ', 1), '_[0-9]+$', '') AS contig_name,
                split_part(Column0, ' ', 1) AS contig_name_6tr,
                Column1 AS aa_seq,
                length(Column1) AS aa_length
            FROM read_csv('/dev/stdin')
            ORDER BY
                contig_name,
                contig_name_6tr
        ) TO '{output.fn_seqs_parq}' (FORMAT 'PARQUET', CODEC 'ZSTD')
EOF
        echo "SQL command"
        cat "$AA_ORFS_SEQS_SQL"
        cat {input.fns_seqs} | seqkit fx2tab | duckdb :memory: -f "$AA_ORFS_SEQS_SQL"
        rm -f "$AA_ORFS_SEQS_SQL"
        """
