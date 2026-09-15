rule get_counts_col_names:
    input:
        fn_merged_all = lambda w: get_merged_all_for_meta(w.batch_meta)
    output:
        fmt_counts_col_names,
    params:
        search = lambda w: get_meta_tab_val(w.batch_meta, 's'),
    shell:
        """
        duckdb \
            --noheader \
            -c ".mode list" \
            -c "
                SELECT column_name 
                FROM (DESCRIBE SELECT * FROM '{input}')
                WHERE column_name LIKE '%{params.search}%'
            " \
            > {output:q}
        """

# duckdb --noheader -c ".mode list" -c " SELECT column_name FROM (DESCRIBE SELECT * FROM '');"

rule get_metadata:
    input:
        meta = lambda w: get_meta_tab_val(w.batch_meta, 'm'),
        norm = lambda w: get_meta_tab_val(w.batch_meta, 'n'),
        fns_count = fmt_counts_col_names,
    output:
        fmt_metadata_batch,
    params:
        script = f"{config['dir_scripts']}/parse_metadata_table.py",
        extra = lambda w: get_meta_tab_val(w.batch_meta, 'e'),
    shell:
        """
        ARG_E="-e {params.extra}"
        if [[ -z "{params.extra}" ]]; then
            ARG_E=""
        fi
        python {params.script:q} \
            -m {input.meta:q} \
            -n {input.norm:q} \
            -c {input.fns_count:q} \
            -k {wildcards.batch_meta} \
            -o {output:q} \
            $ARG_E 
        """

rule merge_metadata:
    input:
        expand(fmt_metadata_batch, batch_meta=BATCHES_META),
    output:
        fn_metadata_merge,
    run:
        df_out = [pd.read_csv(fn) for fn in input]
        df_out = pd.concat(df_out, axis=0, ignore_index=True)
        df_out.to_csv(output[0], index=False)
