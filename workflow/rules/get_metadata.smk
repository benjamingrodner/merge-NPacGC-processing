
rule get_metadata:
    input:
        m = lambda w: get_input_table_value(w.batch_meta, file_tab_meta, 'fn_metadata'),
        n = lambda w: get_input_table_value(w.batch_meta, file_tab_meta, 'fn_norm_factors'),
        c = lambda w: get_input_table_value(w.batch_meta, file_tab_meta, 'fn_kallisto_names'),
        e = lambda w: get_input_table_value(w.batch_meta, file_tab_meta, 'fn_extra_metadata'),
    output:
        