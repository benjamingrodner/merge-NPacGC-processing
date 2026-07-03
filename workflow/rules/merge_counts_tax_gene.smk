rule merge_counts_tax_gene:
    input:
        fn_merged_counts = fmt_merged_counts,
        fn_tax_lin = fmt_tax_lin,
        fn_genes = lambda w: get_fn_genes(w.batch),
    output:
        fn_merged_all = fmt_merged_all
    threads: 21
    resources:
        mem="350G"
    run:
        batch = wildcards.batch
        join_key_counts = config['keys_kallisto']['join']
        join_key_tax = config['keys_diamond']['join']
        join_key_genes = get_input_table_value(
            batch, input_table, 'colname_contig_ko'
        )
        # Set up ibis
        con = ibis.duckdb.connect(memory_limit=resources.mem, threads=threads)
        ibis.set_backend(con)
        # load files
        t_counts = ibis.read_parquet(input.fn_merged_counts)
        t_tax = ibis.read_parquet(input.fn_tax_lin)
        # Merge counts and tax lin
        t_merged = t_counts.join(
            t_tax, 
            predicates=t_counts[join_key_counts] == t_tax[join_key_tax],
            how="outer"
        )
        # Custom load and filter for NS kofam files concatenated with column names intact
        types_ko = {}
        if re.search(r'^G\dNS', batch) is not None:
            print(f"\n\n{ibis.__version__}\n\n")
            types_eval = {"kofam_eval": "string"}
            t_genes = ibis.read_csv(input.fn_genes, types=types_eval)
        else:
            t_genes = ibis.read_csv(input.fn_genes) 
        if re.search(r'^G\dNS', batch) is not None:
            t_genes = t_genes.filter(
                t_genes['kofam_eval'] != 'E-value'
            ).mutate(
                kofam_eval = t_genes['kofam_eval'].try_cast('double')
            )
        # t_counts = t_counts.limit(10000)
        # t_counts = t_counts.limit(10000)
        # t_genes = t_genes.limit(10000)
        # Merge counts and tax lin files
        # Set up taxid enum type
        def create_enum_type(enum_type, vals):
            enum_vals = "', '".join([str(t) for t in vals_sub])
            enum_vals = "('" + enum_vals + "')"
            # with open(f'enumvals_{enum_type}.txt','w') as f:
            #     f.write(enum_vals)
            con.raw_sql(
                f"CREATE TYPE {enum_type} AS ENUM {enum_vals}"
            )
        # Set up to make kofam columns enum categories
        sql_cast = f"SELECT {join_key_genes}, "
        dict_col_caseexpr = {}
        cols_toenum = get_input_table_value(
            batch, input_table, 'colnames_toenum_ko'
        ).split(',')
        for col in cols_toenum:
            # Get unique values
            kos_unique = t_genes.select(col).distinct().execute()
            kos_unique = kos_unique[col].tolist()
            # remove quotes from the values
            vals_sub = []
            case_branches = []
            for v in kos_unique:
                v_ = re.sub("'","(prime)",v)
                v_ = re.sub('"','(2prime)',v_)
                vals_sub.append(v_)
                if v_ != v:
                    case_branches.append((t_genes[col] == v, v_))
            if case_branches:
                dict_col_caseexpr[col] = ibis.cases(
                    *case_branches, else_=t_genes[col]
                )
            # Make enum type
            enum_type = f'{col}_type'
            create_enum_type(enum_type, vals_sub)
            sql_cast += f"{col}::{enum_type} as {col}, "
        # substitute quotation marks in columns
        if dict_col_caseexpr:
            t_genes = t_genes.mutate(**dict_col_caseexpr)
        # Subset to only the columns wanted
        cols_other = get_input_table_value(
            batch, input_table, 'colnames_other_ko'
        ).split(',')
        for col in cols_other:
            sql_cast += f"{col}, "
        # Cast columns to enum and subset
        alias_ibis = 'lifetheuniverseandeverything'
        sql_cast += f"FROM {alias_ibis}"
        t_genes_cat = t_genes.alias(alias_ibis).sql(sql_cast)
        # remove frame selection from tail end of string
        # also keep the frame selected name
        t_genes_cat = t_genes_cat.mutate(
            **{
                'contig_name_6tr': t_genes_cat[join_key_genes],
                join_key_genes: t_genes_cat[join_key_genes].re_replace(r'_\d+$', ''),
            }
        )
        # Merge genes table
        t_merged = t_merged.join(
            t_genes_cat,
            predicates=t_counts[join_key_counts] == t_genes_cat[join_key_genes],
            how="outer"
        )
        # Select columns and save
        columns_select = list(t_counts.columns) # use the counts join key as the contig column
        for t, k in zip([t_tax, t_genes_cat],[join_key_tax, join_key_genes]):
            for col in t.columns:
                if col != k: # remove only the join key from the taxa
                    columns_select.append(col)
        dict_renames = {} # Make NS naming like PA naming
        for c, rn in zip(['kofam','kofam_eval'],['KO','E-value']):
            if c in t_merged.columns:
                dict_renames[rn] = c
        selected = t_merged.select(columns_select).rename(**dict_renames)
        print('Merging counts, taxa lineages, and gene annnotations.')
        selected.to_parquet(output.fn_merged_all)


rule check_stuff:
    input:
        fn_merged_counts = fmt_merged_counts,
        fn_merged_all = fmt_merged_all
    output:
        testfile = fmt_testfile
    threads: 1
    resources:
        mem="1G"
    run:
        # Set up ibis
        con = ibis.duckdb.connect(memory_limit=resources.mem, threads=threads)
        ibis.set_backend(con)
        # load files
        t = ibis.read_parquet(input.fn_merged_all)
        with open(output.testfile, 'w') as f:
            # check a specific target id
            tid = t.select(config['keys_kallisto']['join']).head(1).to_pandas().iloc[0, 0]
            dfsub = t.filter(t[config['keys_kallisto']['join']] == tid).to_pandas()
            for i, row in dfsub.iterrows():
                # f.write(f'\n{i}')
                for k, v in row.items():
                    f.write(f'{k},{v},\n')
            # # check values
            # df = t.limit(3).to_pandas()
            # for i in range(df.shape[0]):
            #     for k, v in df.iloc[i,:].squeeze().items():
            #         f.write(f'{v},{k},\n')
        

