# =============================================================================
# Functions
# =============================================================================

def get_tar_names(fn_tar, re_tar, batch):
    # Make tar name dir
    out_dir = config['path_tarnames'] + '/tar_names'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    # Make filename
    fn_tar_names = os.path.split(fn_tar)[1]
    fn_tar_names = re.sub('.tar.gz','-tar_names.txt', fn_tar_names)
    fn_tar_names_full = f'{out_dir}/{batch}-{fn_tar_names}'
    # Check if the file already exists
    print(fn_tar_names_full)
    if not os.path.exists(fn_tar_names_full):
        print('Getting tar sub-filenames')
        tar_names = []
        # Read the tar names
        with tarfile.open(fn_tar, "r:gz") as tar:
            for name in tar:
                n = name.name
                if bool(re.search(re_tar, n)):
                    tar_names.append(n)
        # Write file
        with open(fn_tar_names_full, 'w') as f:
            for line in tar_names:
                f.write(f"{line}\n")
    # If it already exists, read the names
    else:
        with open(fn_tar_names_full, 'r') as f:
            tar_names = f.read().splitlines()
    return tar_names

def get_glob_fns(dir_glob, fn_glob):
    if dir_glob:
        return(glob.glob(f'{dir_glob}/{fn_glob}'))
    else:
        return(glob.glob(fn_glob))



def get_fns(fmt, col_dir='search_output_dir', col_fn_metat='fn_metat', col_fn_targets='fn_targets', col_check=''):
    fns = []
    for index, row in input_table.iterrows():
        fn = fmt.format(
            search_output_dir=row[col_dir],
            fn_metat=row[col_fn_metat], 
            fn_targets=row[col_fn_targets],
        )
        if not col_check:
            # Build and append filename
            fns.append(fn)
        else:
            if row[col_check]:
                fns.append(fn)

    return fns

def get_fns_sample(fmt, col_check=''):
    fns = []
    for index, row in input_table.iterrows():
        # Column to check if appending row filenames
        appnd = False
        if not col_check:
            # Build and append filename
            appnd = True
        else:
            if row[col_check]:
                appnd = True
        if appnd:
            # Use the glob expression to get a list of kallisto filenames
            fns_kallisto = get_glob_fns(
                row.dir_kallisto, row.glob_kallisto
            )
            for fnk_ in fns_kallisto:
                # Check if the file is a tarball
                if row.re_kallisto_tar:
                    # Extract the filenames in the tarball 
                    tar_names = get_tar_names(
                        fnk_, 
                        row.re_kallisto_tar, 
                        row.batch
                    )
                    for fnkt_ in tar_names:  
                        fnkt = os.path.split(fnkt_)[1]
                        fnk = os.path.split(fnk_)[1]
                        fnk_merge = f'{fnk}.{fnkt}'
                        # Build and append output filename
                        fns.append(fmt.format(
                            batch=row.batch, 
                            fn_kallisto=fnk_merge
                        ))
                else: 
                    fnk = os.path.split(fnk_)[1]
                    # Build and append filename
                    fns.append(fmt.format(
                        batch=row.batch,
                        fn_kallisto=fnk
                    ))
    return fns

def get_fns_row_sample(
    fmt, 
    batch, input_table, 
    col_dir='search_output_dir', 
    col_fn_metat='fn_metat', 
    col_fn_targets='fn_targets', 
    col_check=''
):
    fns = []
    row = get_input_table_row(batch, input_table)
    # Column to check if appending row filenames
    appnd = False
    if not col_check:
        # Build and append filename
        appnd = True
    else:
        if row[col_check]:
            appnd = True
    if appnd:
        # Use the glob expression to get a list of kallisto filenames
        fns_kallisto = get_glob_fns(
            row.dir_kallisto, row.glob_kallisto
        )
        for fnk_ in fns_kallisto:
            # Check if the file is a tarball
            if row.re_kallisto_tar:
                # Extract the filenames in the tarball 
                tar_names = get_tar_names(
                    fnk_, 
                    row.re_kallisto_tar, 
                    row.search_output_dir,
                    row.fn_metat
                )
                for fnkt_ in tar_names:  
                    fnkt = os.path.split(fnkt_)[1]
                    fnk = os.path.split(fnk_)[1]
                    fnk_merge = f'{fnk}.{fnkt}'
                    # Build and append output filename
                    fns.append(fmt.format(
                        search_output_dir=row.search_output_dir,
                        fn_metat=row.fn_metat, 
                        fn_targets=row.fn_targets,
                        fn_kallisto=fnk_merge
                    ))
            else: 
                fnk = os.path.split(fnk_)[1]
                # Build and append filename
                fns.append(fmt.format(
                    search_output_dir=row.search_output_dir,
                    fn_metat=row.fn_metat, 
                    fn_targets=row.fn_targets,
                    fn_kallisto=fnk
                ))

    return fns

def get_fns_target(fmt, col_dir='search_output_dir', col_fn_targets='fn_targets', col_check=''):
    fns = []
    for index, row in input_table.iterrows():
        fn = fmt.format(
            search_output_dir=row[col_dir],
            fn_targets=row[col_fn_targets],
        )
        if not col_check:
            # Build and append filename
            fns.append(fn)
        else:
            if row[col_check]:
                fns.append(fn)
    return set(fns)

def check_contigs_from_other(batch, input_table):
    contigs_from_other = get_input_table_value(
        batch, input_table, 'get_contigs_from_other_hmmsearch'
    )
    return  contigs_from_other if contigs_from_other else batch

# To avoid circular definitino in check contigs from other
def get_input_table_value(batch, input_table, column):
    return input_table.loc[(input_table.batch == batch), column].values[0]

#  function to get value checks for batch substitution
def get_input_table_value_subbatch(batch, input_table, column):
    batch = check_contigs_from_other(batch, input_table)
    return input_table.loc[(input_table.batch == batch), column].values[0]

def get_input_table_row(batch, input_table):
    return input_table.loc[(input_table.batch == batch), :].squeeze()


def get_fn_full_subbatch(batch, input_table, col_dir, col_fn):
    batch = check_contigs_from_other(batch, input_table)
    path = get_input_table_value(batch, input_table, col_dir)
    fn = get_input_table_value(batch, input_table, col_fn)
    if path:
        return path + f'/{fn}'
    else:
        return fn

def get_fn_full(batch, input_table, col_dir, col_fn):
    path = get_input_table_value(batch, input_table, col_dir)
    fn = get_input_table_value(batch, input_table, col_fn)
    if path:
        return path + f'/{fn}'
    else:
        return fn


def get_fn_kallisto_full(batch, fn_kallisto, input_table):
    row = get_input_table_row(batch, input_table)
    # Split wildcard to get kallisto tarball filename
    if row.re_kallisto_tar:
        regex = fnmatch.translate(row.glob_kallisto)
        regex = regex.rstrip('\Z')
        fn_kallisto = re.match(regex, fn_kallisto)[0]
    # Assemble full filename
    fns_glob = get_glob_fns(
        row.dir_kallisto, row.glob_kallisto
    )
    for fnk_full in fns_glob:
        if fn_kallisto in fnk_full:
            return fnk_full

def get_fn_kallisto_tar(batch, fn_kallisto, input_table):
    row = get_input_table_row(
        batch, input_table,
    )
    # Split wildcard to get kallisto tarball filename
    if row.re_kallisto_tar:
        regex = fnmatch.translate(row.glob_kallisto + '.')
        regex = regex.rstrip('\Z')
        fn_kallisto_tar = re.split(regex, fn_kallisto)[1]
        return fn_kallisto_tar
    # Assemble full filename
    else:
        return ''

def get_fns_counts(batch):
    row = get_input_table_row(batch, input_table)
    fnc = f"{row['dir_kallisto']}/{row['glob_kallisto']}"
    return glob.glob(fnc)

def get_fn_tax(batch):
    row = get_input_table_row(batch, input_table)
    return f"{row['dir_diamond']}/{row['fn_diamond']}"

def get_fn_genes(batch):
    row = get_input_table_row(batch, input_table)
    return f"{row['path_kofam']}/{row['fn_kofam']}"

# Parse the kallisto sample name (the sub tar name in tarball)
def get_meta(sn_type):
    meta_fn = input_table.loc[
        input_table['sn_type_parse_kallisto'] == sn_type, 
        'fn_sample_metadata'
    ].values[0]
    return pd.read_csv(meta_fn)

def get_size_lat_depth_rep_timep(meta_sample, rnd_lat, skip=[]):
    vals = []
    if not 'size' in skip:
        size = str(float(meta_sample['Filter'].values[0]))
        size += 'um'
        vals.append(size)
    if not 'lat' in skip:
        lat = str(round(float(meta_sample['Latitude'].values[0]), rnd_lat))
        lat += 'deg'
        vals.append(lat)
    if not 'depth' in skip:
        depth = str(float(meta_sample['Depth'].values[0]))
        depth += 'm'
        vals.append(depth)
    if not 'rep' in skip:
        rep = meta_sample['Replicate'].values[0]
        vals.append(rep)
    if not 'timep' in skip:
        timep = meta_sample['Datetime'].values[0]
        timep = re.sub('/','_', timep)
        timep = re.sub(r'\s','-', timep)
        vals.append(timep)
    return vals

def parse_fn_kallisto_sn(fn='', sn_type='', get_columns=False, rnd_lat=2):
    if not get_columns:
        ass, sample, lat, ammend, timep, depth, size, rep = [''] * 8
        if sn_type == 'G1NS':
            splt = fn.split('.')
            ass, sm_sz, rp = splt[:3]
            sample, sz = sm_sz.split('_',1)
            # size = str(float(re.sub('_','.',sz)))
            alias2 = sm_sz + rp
            meta = get_meta(sn_type)
            meta_sample = meta.loc[meta['Alias2'] == alias2, :]
            size, lat, depth, rep, timep = get_size_lat_depth_rep_timep(meta_sample, rnd_lat)
        elif sn_type == 'G2NS':
            ass, sample, dp, sz, rp, _ = fn.split('.')
            meta = get_meta(sn_type)
            alias2 = f'{sample}.{dp}.{sz}.{rp}'
            meta_sample = meta.loc[meta['Alias2'] == alias2, :]
            size, lat, depth, rep, timep = get_size_lat_depth_rep_timep(meta_sample, rnd_lat)
        elif sn_type == 'G3NS':
            meta = get_meta(sn_type)
            sid = os.path.splitext(fn)[0]
            meta_sample = meta.loc[meta['SampleID'] == sid, :]
            size, lat, depth, rep, timep = get_size_lat_depth_rep_timep(meta_sample, rnd_lat)
            ass = re.match(r'.+NS', fn)[0]
            sample = re.search(r'UW\d+_\d', fn)[0]
        elif sn_type == 'G5':
            ass, sample, ammend, timep, rep, _ = fn.split('.')
        elif sn_type == 'D1':
            ass, sm_rep_tp, _, _ = fn.split('.')
            sample, rep, timep = sm_rep_tp.split('_')
        elif sn_type == 'G1PA':
            meta = get_meta(sn_type)
            ass, fn_ = fn.split('.', 1)
            sid = re.match(r'.+(?=\.abundance)', fn_)[0]
            sid = re.sub(r'\.','_',sid)
            meta_sample = meta.loc[meta['SampleID'] == sid, :]
            size, lat, depth, rep, timep = get_size_lat_depth_rep_timep(meta_sample, rnd_lat)
            sample, _ = fn_.split('_', 1)
        elif sn_type == 'G2PA':
            _, ass, sample, dp, sz, rp, _, _ = fn.split('.')
            meta = get_meta(sn_type)
            sid = f"{sample}.{dp}.{sz}.{rp}"
            meta_sample = meta.loc[meta['SampleID'] == sid, :]
            size, lat, rep = get_size_lat_depth_rep_timep(
                meta_sample, 
                rnd_lat, 
                skip=['depth','timep']
            )
            depth = str(float(dp[:-1])) + 'm'
            # size = re.sub('_','.',sz)
        elif sn_type == 'G3PA.UW':
            meta = get_meta(sn_type)
            ass, sample_ = fn.split('.')[:2]
            meta_sample = meta.loc[meta['Alias2'] == sample_, :]
            size, lat, depth, rep, timep = get_size_lat_depth_rep_timep(meta_sample, rnd_lat)
            sample_list = str(meta_sample['Alias1'].values[0]).split(' ')
            sample = sample_list[0]
            if '#' in sample_list[1]:
                sample += sample_list[1]
        elif sn_type == 'G3PA.diel':
            ass1, ass2, sample, rp, _, _, _, _ = fn.split('.')
            ass = f'{ass1}.{ass2}'
            meta = get_meta(sn_type)
            sid = f"{ass}.{sample}.{rp}"
            meta_sample = meta.loc[meta['SampleID'] == sid, :]
            size, lat, depth, rep, timep = get_size_lat_depth_rep_timep(meta_sample, rnd_lat)

        elif sn_type == 'G3PA.PM':
            ass = sn_type
            sample = re.search(r'UW\d+_\d', fn)[0]
            meta = get_meta(sn_type)
            sid = re.match(r'.+(?=\.unstranded)', fn)[0]
            meta_sample = meta.loc[meta['SampleID'] == sid, :]
            size, lat, depth, rep, timep = get_size_lat_depth_rep_timep(meta_sample, rnd_lat)
        else:
            raise ValueError(
                f"""
                Sample name parse type {sn_type} not configured or not provided 
                (sn_type_parse_kallisto column in file table)
                """
            )
        return [ass, sample, lat, ammend, timep, depth, size, rep]
    else:
        return ['assembly', 'sample', 'latitude','ammendment', 'timepoint', 'depth', 'size', 'rep']

def get_fn_list_contigs_diamond(batch, input_table, fn_6tr, fn_std):
    type_contig_name = get_input_table_value(
            batch, input_table, 'type_diamond_contig_name'
        )
    if type_contig_name == '6tr':
        return fn_6tr
    else: 
        return fn_std

def get_fns_batch(fmt):
    return [fmt.format(batch=batch) for batch in config['batches_to_run']]

def calculate_cores(frac=1):
    return max(1, int(workflow.cores * frac))

def calculate_mem(frac=1):
    global_mem = workflow.resource_settings.get('mem_gb')
    if global_mem is not None:
        return max(1024, int(global_mem * frac)) # Use at least 1 GB
    else:
        return 4096 # Default to 4 GB if --resources was not specified
# def get_fn_metat_full(fn_metat, fn_targets, input_table):
#     path = get_input_table_value(fn_metat, fn_targets, input_table, 'path_metat')
#     if path:
#         return path + f'/{fn_metat}'
#     else:
#         return fn_metat


# def get_fn_targets_full(fn_metat, fn_targets, input_table):
#     path = get_input_table_value(fn_metat, fn_targets, input_table, 'path_targets')
#     if path:
#         return path + f'/{fn_targets}'
#     else:
#         return fn_targets
# Get sequence filenames
def expand_fastanames(fmt, dict_batch_fastaname_fnseqs, batch):
    fns = []
    for f, _ in dict_batch_fastaname_fnseqs[batch].items():
        fn = fmt.format(batch="{batch}", ko="{ko}",fastaname=f)
        fns.append(fn)
    return fns


def get_fnseqs(batch, input_table):
    glob_seqs = get_input_table_value(batch, input_table, 'glob_sequences')
    return glob.glob(glob_seqs)
    

