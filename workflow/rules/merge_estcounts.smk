
rule merge_estcounts:
    input:
        fns_counts = lambda w: get_fns_counts(w.batch)
    output:
        fn_merged_counts = temp(fmt_merged_counts) 
        # fn_all_taxon_estcounts = fmt_all_taxon_estcounts
    threads: 42
    resources:
        mem="700G"
    run:
        # Set up ibis
        con = ibis.duckdb.connect(memory_limit=resources.mem, threads=threads)
        ibis.set_backend(con)
        print(f'Set up duckdb with mem={resources.mem} threads={threads}')
        # Filename info
        input_table = pd.read_csv(
            config['input_table'],
            keep_default_na=False
        )
        batch = wildcards.batch
        # Get a list of files and whether or not the counts are merged yet
        mergedyn = get_input_table_value(
            batch, input_table, 'merged_kallisto_yn'
        )
        dict_fnc_mergedyn = {fnc: mergedyn for fnc in input.fns_counts}
        # Get a list of tarnames
        dict_fnc_tarnames = {}
        for fnc in input.fns_counts:
            if '.tar.gz' in fnc:
                re_tar = get_input_table_value(
                    batch, input_table, 're_kallisto_tar'
                )
                tarnames = get_tar_names(fnc, re_tar, batch)
                dict_fnc_tarnames[fnc] = tarnames


        # Add experiment files that are mapped to the batch
        bexps = config['dict_batch_exps'].get(batch)
        if bexps is not None:
            for batchexp in bexps:
                row = get_input_table_row(batchexp, input_table)
                fnc_glob = f"{row['dir_kallisto']}/{row['glob_kallisto']}"
                for fnc in glob.glob(fnc_glob):
                    dict_fnc_mergedyn[fnc] = row['merged_kallisto_yn']
                    # Check if tarball and add tarnames
                    if '.tar.gz' in fnc:
                        re_tar = get_input_table_value(
                            batchexp, input_table, 're_kallisto_tar'
                        )
                        tarnames = get_tar_names(fnc, re_tar, batchexp)
                        dict_fnc_tarnames[fnc] = tarnames

        # Set up merge
        join_key = config['keys_kallisto']['join']
        counts_key = config['keys_kallisto']['counts']
        t_merged = None
        cols_select = [join_key]
        for fnc, mergedyn in dict_fnc_mergedyn.items():
            # Load file
            tarnames = dict_fnc_tarnames.get(fnc)
            if tarnames is None:
                t = ibis.read_csv(fnc)
                if mergedyn == 'n':
                    # Rename columns with sample name id
                    sn = os.path.split(fnc)[-1]
                    t = t.rename({sn: counts_key})
                    cols_select.append(sn)
                else:
                    cols_select += [c for c in t.columns if c != join_key]
            else: # load from tarball
                print('Loading from tarball...')
                with tarfile.open(fnc, "r:gz") as tar:
                    t = None
                    for tnm in tarnames:
                        print(tnm)
                        f = tar.extractfile(tnm)
                        parse_options = pv.ParseOptions(delimiter="\t")
                        read_options = pv.ReadOptions(block_size=20_000_000)
                        if os.path.splitext(tnm)[1] == '.gz':
                            with gzip.GzipFile(fileobj=f, mode='rb') as f_decomp:
                                reader = pv.open_csv(
                                    f_decomp, 
                                    read_options=read_options, 
                                    parse_options=parse_options
                                )
                                pa_table = reader.read_all()
                                t_tnm = con.create_table(tnm, pa_table)
                        else:
                            reader = pv.open_csv(
                                f, 
                                read_options=read_options, 
                                parse_options=parse_options
                            )
                            pa_table = reader.read_all()
                            t_tnm = con.create_table(tnm, pa_table)


                        if mergedyn == 'n':
                            # Rename columns with sample name id
                            sn = os.path.split(tnm)[-1]
                            t_tnm = t_tnm.rename({sn: counts_key})
                            cols_select.append(sn)
                        else:
                            cols_select += [c for c in t_nm.columns if c != join_key]
                        if t is None:
                            t = t_tnm
                        else:
                            t = t.join(
                                t_tnm,
                                predicates=t[join_key] == t_tnm[join_key],
                                how="outer" # 
                            )
            # Merge tables
            print(f'Merging {fnc}')
            if t_merged is None:
                t_merged = t
            else:
                t_merged = t_merged.join(
                    t,
                    predicates=t_merged[join_key] == t[join_key],
                    how="outer" # 
                )
        selected = t_merged.select(cols_select)
        selected.to_parquet(output.fn_merged_counts)