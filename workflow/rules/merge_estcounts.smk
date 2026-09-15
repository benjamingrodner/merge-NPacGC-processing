
# for files where the sample name is stored in the tmp column
def split_count_file(fn, dir_out, batch):
    header = 'target_id\tlength\teff_length\test_counts\ttpm|name_file\n'
    dict_sn_handle = {}
    
    os.makedirs(dir_out, exist_ok=True)
    dict_sn_fn = {}
    fmt_out = f'{dir_out}/{batch}-kallisto_counts-{{sn}}'
    
    dir_done = f'{dir_out}/done'
    os.makedirs(dir_done, exist_ok=True)
    fn_done = f'{dir_done}/{batch}-split_done.txt'
    
    if not os.path.exists(fn_done):
        print('Splitting file:', batch)
        i = 0
        with gzip.open(fn, 'rt') as f:
            for line in f:
                l = line.strip().split('\t')
                if (l[0] != 'target_id') and (len(l[0]) > 0):
                    try:
                        sn = l[4].split('|')[1]
                    except:
                        raise IndexError(f'Line = {line}\nl = {l}')
                    fn_out = fmt_out.format(sn=sn)
                    handle = dict_sn_handle.get(sn)
                    if handle is None:
                        handle = open(fn_out, 'w')
                        handle.write(header)
                        dict_sn_handle[sn] = handle
                    else:
                        handle = dict_sn_handle[sn]
                    handle.write(line)
                i += 1
                if i % 1e6 == 0:
                    print(f'Lines read: {i}', end='\r')

        for sn, handle in dict_sn_handle.items():
            handle.close()
        # done file
        with open(fn_done, "w") as f:
            f.write('done')
    else:
        print('File already split:', batch)

    return fmt_out.format(sn='*'), fn_done


rule merge_estcounts:
    input:
        fns_counts = lambda w: get_fns_counts(w.batch)
    output:
        fn_merged_counts = fmt_merged_counts
        # fn_all_taxon_estcounts = fmt_all_taxon_estcounts
    threads: 32
    resources:
        mem="200G"
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
            batch, input_table, 'merged_kallisto_ynm'
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
        split_counts_fns = []
        if bexps is not None:
            for batchexp in bexps:
                row = get_input_table_row(batchexp, input_table)
                fnc_glob = f"{row['dir_kallisto']}/{row['glob_kallisto']}"
                for fnc in glob.glob(fnc_glob):
                    ynm = row['merged_kallisto_ynm']
                    if ynm == 'm': # for files that need splitting
                        fncg, fn_split_done = split_count_file(
                            fnc, config['split_count_dir'], batchexp
                        )
                        split_counts_fns.append(fn_split_done)
                        for fncs in glob.glob(fncg):
                            split_counts_fns.append(fncs)
                            dict_fnc_mergedyn[fncs] = 'n'
                    else:
                        dict_fnc_mergedyn[fnc] = ynm
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
        
        # Remove temp kallisto split
        for fn in split_counts_fns:
            os.remove(fn)