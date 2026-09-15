def filter_6tr_rows(t, group_col, order_cols, ascending):
    order_exprs = [
        t[col].asc() if asc else t[col].desc()
        for col, asc in zip(order_cols, ascending)
    ]

    w = ibis.window(
        group_by=t[group_col],
        order_by=order_exprs,
    )
    return (
        t
        .mutate(rn=ibis.row_number().over(w))
        .filter(lambda t: t.rn == 0)
        .drop("rn")
    )


rule get_tax_lineage:
    input:
        fn_tax = lambda w: get_fn_tax(w.batch)
    output:
        fn_tax_lin = fmt_tax_lin
    threads: 8
    resources:
        mem="50G"
    run:
        # Set up ibis
        con = ibis.duckdb.connect(memory_limit=resources.mem, threads=threads)
        ibis.set_backend(con)
        # Load file and give column names
        table_name = 'table_diamond'
        t = ibis.read_csv(
            input.fn_tax, 
            header=False, 
            names=config['keys_diamond']['cols'],
            table_name=table_name
            )
        print('COL T',t.columns)
        # t = t.limit(100)
        # Get unique taxids
        print('Getting unique taxa')
        # taxids = t['taxid'].cast('string')
        taxid_key = config['keys_diamond']['taxid']
        taxa_unique = t.select(taxid_key).distinct().execute()
        taxa_unique = taxa_unique[taxid_key].tolist()
        print(f'There are {len(taxa_unique)} unique taxa')
        print('Getting dict of taxa lineages')
        # Make a dict mapping taxids to lineage
        ncbi = NCBITaxa()
        cols_lin = config['cols_lin']
        dict_col_tax_nm = defaultdict(dict)
        dict_col_nms = defaultdict(list)

        for tax in taxa_unique:
            if tax != 0:
                try:
                    lin_full = ncbi.get_lineage(tax)
                except:
                    lin_full = None
                if lin_full is not None:
                    dict_t_n = ncbi.get_taxid_translator(lin_full)
                    ranks = {r:tid for tid,r in ncbi.get_rank(lin_full).items()}
                    lin_t = [ranks.get(c, '') for c in cols_lin]
                    lin_nm = [dict_t_n.get(tid, '') for tid in lin_t]
                    for col, nm in zip(cols_lin, lin_nm):
                        dict_col_tax_nm[col][tax] = nm
                        dict_col_nms[col].append(nm)
                # else:
                #     for col in cols_lin:
                #         dict_col_tax_nm[col][tax] = ''
            # else:
            #     for col in cols_lin:
            #         dict_col_tax_nm[col][tax] = ''

        # Create new columns with mapping logic
        dict_col_expr = {}
        for col, dict_tax_nm in dict_col_tax_nm.items():
            case_branches = []
            for tax, nm in dict_tax_nm.items():
                case_branches.append((t[taxid_key] == tax, nm))
            cases_expr = ibis.cases(*case_branches, else_='').cast("string")
            dict_col_expr[col] = cases_expr
        t_mapped = t.mutate(**dict_col_expr)
        print('COL EXPR',t_mapped.columns)
        # Set up taxid enum type
        def create_enum_type(enum_type, vals):
            enum_vals = "', '".join([str(v) for v in vals])
            enum_vals = "('" + enum_vals + "')"
            # with open(f'enumvals_{enum_type}.txt','w') as f:
            #     f.write(enum_vals)
            con.raw_sql(
                f"CREATE TYPE {enum_type} AS ENUM {enum_vals}"
            )
        enum_type = f'{taxid_key}_type'
        create_enum_type(enum_type, taxa_unique)
        # Set up sql cast 
        join_key = config['keys_diamond']['join']
        sql_cast = f"SELECT {join_key}, {config['keys_diamond']['col_eval']}, "
        sql_cast += f"{taxid_key}::{enum_type} as {taxid_key}, "
        dict_col_caseexpr = {}
        # Set up tax name enum types and add to cast statement
        for col, nms in dict_col_nms.items():
            enum_type = f'{col}_type'
            nms = list(set(nms))
            # remove quotes from the values
            vals_sub = []
            case_branches = []
            for v in nms:
                v_ = re.sub("'","(prime)",v)
                v_ = re.sub('"','(2prime)',v_)
                vals_sub.append(v_)
                if v_ != v:
                    case_branches.append((t_mapped[col] == v, v_))
            if case_branches:
                dict_col_caseexpr[col] = ibis.cases(
                    *case_branches, else_=t_mapped[col]
                )
            # print(f'Creating enum type {enum_type}')
            create_enum_type(enum_type, vals_sub)
            sql_cast += f"NULLIF({col}, '')::{enum_type} as {col}, "
        # substitute quotation marks in columns
        if dict_col_caseexpr:
            t_mapped = t_mapped.mutate(**dict_col_caseexpr)
        print('COL CASE',t_mapped.columns)
        # TODO Set up a column with the full lineage as a string
        # Mutate columns to enum categories
        alias_ibis = 'lifetheuniverseandeverything'
        sql_cast += f"FROM {alias_ibis}"
        colstrcast = [taxid_key] + cols_lin
        dict_strcast = {col: 'string' for col in colstrcast}
        t_mapped = t_mapped.cast(dict_strcast)
        t_cat = t_mapped.alias(alias_ibis).sql(sql_cast)
        print('COL ALIAS',t_cat.columns)
        # t_mapped = t_mapped.mutate(
        #     **{
        #         col: t_mapped[col].cast(enum_type) 
        #         for col, enum_type in dict_col_enumtype.items()
        #     }
        # )
        # remove frame selection from tail end of string
        t_cat = t_cat.mutate(**{
            config['keys_diamond']['col_6tr']: t_cat[join_key],
            join_key: t_cat[join_key].re_replace(r'_\d+$', '')
        })
        print('COL 6TR',t_cat.columns)
        # remove duplicates for join key, selecting frame with best KO hit
        t_cat = filter_6tr_rows(
            t_cat, join_key, [config['keys_diamond']['col_eval']], [True]
        )
        # Check for empty string in 
        # Write to parquet
        print('Writing to file')
        t_cat.to_parquet(output.fn_tax_lin)

