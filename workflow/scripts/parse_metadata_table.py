# skeleton produced using gemini 3.5 flash in free tier 5/31/26
### prompt
# Please write a python script with click arguments:
# - a filename to a csv table with metadata
# - a filename to another csv table with normalization factors
# - a string that is a key to decide which function to use to parse the table
# - an output filename

# Steps:
# - load the csvs to dataframes using pandas
# - use the string to select the correct function
# - execute the function on the dataframes which produces a  new dataframe
# - write the new dataframe to the output filename
###

# Edits Ben Grodner

# Refactor funcs gemini 3.5 flash free tier

import re
import sys
import click
import pandas as pd


# -----------------------------------------------------------------------------
# Parsing Functions
# -----------------------------------------------------------------------------

import re
import pandas as pd
from typing import Tuple, Optional

def _parse_station(
    df_meta: pd.DataFrame, 
    df_norm: pd.DataFrame, 
    regex: str, 
    col_norm: str,
    key_formatter: callable,
    col_meta = 'SampleID',
    col_meta_merge = 'norm_key',
    ignore_list = []
) -> pd.DataFrame:
    """Core driver function that handles extraction, validation, and merging."""
    df_meta = df_meta.copy()
    
    # Extract common regex groups: S(\d+), (C\d+), (.+um), (\w)
    
    keys_counts = []
    keys_norm = []
    
    for sample_id in df_meta[col_meta]:
        match = re.search(regex, sample_id)
        if not match:
            raise ValueError(f"SampleID '{sample_id}' does not match regex '{regex}' pattern.")
                    
        # Delegate format variations to the specific type formatter
        k, kn = key_formatter(match)
        
        keys_counts.append(k)
        keys_norm.append(kn)
    
    df_meta[COL_KALLISTO_MAP] = keys_counts
    df_meta[col_meta_merge] = keys_norm
    
    # Validation
    _validate_merge(df_meta, df_norm, col_meta_merge, col_norm, ignore_list)
        
    return df_meta.merge(df_norm, how='left', left_on=col_meta_merge, right_on=col_norm)

def _validate_merge(df_meta: pd.DataFrame, df_norm: pd.DataFrame, mcol: str, ncol: str, ignore_list=[]):
    # --- Validation & Merge (Same as before) ---
    meta_ids = set(df_meta[mcol].dropna())
    norm_ids = set(df_norm[ncol].dropna())
    
    meta_miss = list(meta_ids - norm_ids)
    norm_miss = list(norm_ids - meta_ids)
    if ignore_list:
        meta_miss = [k for k in meta_miss if not any([i in k for i in ignore_list])]
        norm_miss = [k for k in norm_miss if not any([i in k for i in ignore_list])]

    if meta_miss or norm_miss:
        raise ValueError(
            f"\nMismatch between key in metadata and merge column:\n"
            f"Missing in merge column (from key {mcol}): {meta_miss}\n"
            f"Missing in key (from merge column {ncol}): {norm_miss}"
        )
    return


# --- Public API Functions ---

def parse_g1ns(df_meta: pd.DataFrame, df_norm: pd.DataFrame) -> pd.DataFrame:
    """Parses G1NS metadata and normalizes."""
    def formatter(match):
        s, c, f, r = match.groups()
        s = s.zfill(2)

        if s == '10': c = ''
        f_underscore = f.replace('.', '_')
        k = f'G1NS.S{s}{c}_{f_underscore}.{r}.tsv'
        kn = f'S{s}{c}_{f_underscore}{r}'
        return k, kn

    prefix = 'g1-st-am-ns'
    regex = rf'(?<={prefix}\.)S(\d+)\.(C\d+)\.(.+um)\.(\w)'

    return _parse_station(df_meta, df_norm, regex, 'sample_name_short', formatter)


def parse_g1pa(df_meta: pd.DataFrame, df_norm: pd.DataFrame) -> pd.DataFrame:
    """Parses G1PA metadata and normalizes."""
    def formatter(match):
        s, c, f, r = match.groups()
        s = s.zfill(2)
        if s == '10': c = ''
        k = f'G1PA.S{s}{c}_{f}{r}.abundance.tsv'
        
        f_clean = f.replace('.', '')
        kn = f'S{s}{c}_{f_clean}{r}'
        
        # Handle the custom S06C1 exception
        if k == 'G1PA.S06C1_3umB.abundance.tsv':
            k = 'G1PA.S06C1_3umC.abundance.tsv'
            kn = 'S06C1_3umC'
        return k, kn

    prefix = 'g1-st-am-pa'
    regex = rf'(?<={prefix}\.)S(\d+)\.(C\d+)\.(.+um)\.(\w)'
    
    return _parse_station(df_meta, df_norm, regex, 'sample_name_short', formatter)


def parse_g2ns_st(df_meta: pd.DataFrame, df_norm: pd.DataFrame) -> pd.DataFrame:
    """Parses G2NS station metadata and normalizes."""
    def formatter(match):
        s, c, f, r = match.groups()
        s = s.zfill(2)
        f_underscore = f.replace('.', '_')
        k = f'G2NS.S{s}{c}.15m.{f_underscore}.{r}.tsv'
        
        kn = f'G2NS.S{s}{c}.15m.{f_underscore}.{r}'
        return k, kn

    prefix = 'g2-st-am-ns'
    regex = rf'(?<={prefix}\.)S(\d+)\.(C\d+)\.(.+um)\.(\w)'

    return _parse_station(df_meta, df_norm, regex, 'sample_name', formatter)


def parse_g2ns_dcm(df_meta: pd.DataFrame, df_norm: pd.DataFrame) -> pd.DataFrame:
    """Parses G2NS DCM metadata and dynamically maps depth from df_meta."""
    df_meta = df_meta.copy()

    mcol = 'SampleID'
    mcol_norm = 'extractionID'
    ncol = 'Sample_ID'

    # Target regex pattern
    regex = r'(?<=g2-ctd-ns\.)S(\d+)\.(C\d+)\.DCM\.(.+um)\.(\w)'
    
    keys_counts = []
    
    # Use zip to iterate over both the SampleID and your Depth column simultaneously
    # Replace 'Depth_m' with the exact name of your depth column
    for sample_id, depth_val in zip(df_meta[mcol], df_meta['Depth.m']):
        match = re.search(regex, sample_id)
        if not match:
            raise ValueError(f"SampleID '{sample_id}' does not match expected DCM pattern.")
            
        s, c, f, r = match.groups()
        s = s.zfill(2)
        
        # Ensure the depth has 'm' attached (e.g., 87 -> '87m' or '87m' -> '87m')
        depth_str = str(depth_val).strip()
        if not depth_str.endswith('m'):
            depth_str = f"{depth_str}m"
        
        # Construct the count key dynamically
        f_underscore = f.replace('.', '_')
        k = f'G2.DCM.NS.S{s}{c}.{depth_str}.{f_underscore}.{r}.abundance.tsv.gz'
        keys_counts.append(k)
        
    df_meta[COL_KALLISTO_MAP] = keys_counts

    # Adjust
    df_norm_adj = df_norm.copy()
    df_norm_adj.loc[df_norm_adj[ncol] == 'BD63', ncol] = 'BD60'
    
    _validate_merge(df_meta, df_norm_adj, mcol_norm, ncol)
        
    return df_meta.merge(df_norm_adj, how='left', left_on=mcol_norm, right_on=ncol)


def parse_g2ns_rexp(df_meta: pd.DataFrame, df_norm: pd.DataFrame) -> pd.DataFrame:
    df_meta = df_meta.copy()
    
    # Target regex pattern
    prefix = 'g2-inc-ns'
    regex = rf'(?<={prefix}\.)S(\d+)\.rexp(\d)\.(\w+)\.(T\d+)\.(.+um)\.(\w)'
    # Specific formatting
    def formatter(match):
        s, e, a, t, f, r = match.groups()
        s = s.zfill(2)
        f_underscore = f.replace('.', '_')
        if (e != "1") & (t == "T72"): t = "T96"

        kn = f'G2.REXP{e}.NS.{a}.{t}.{f_underscore}.{r}'
        k = f'{kn}.abundance.tsv.gz'
        
        return k, kn

    ignore = ['G2.DCM.NS']
    df_norm = df_norm.rename(columns={'SampleID':'SampleID_norm'})
    return _parse_station(
        df_meta, df_norm, regex, 
        col_norm='SampleID_norm', key_formatter=formatter, ignore_list=ignore)

def parse_g2pa_dcm(df_meta: pd.DataFrame, df_norm: pd.DataFrame) -> pd.DataFrame:
    df_meta = df_meta.copy()
    
    mcol = 'extractionID'
    ncol = 'Sample_ID'

    keys_counts = [f"g2-ctd-pa-kallisto_counts-{sid}.abundance.tsv" 
                   for sid in df_meta[mcol]]
        
    df_meta[COL_KALLISTO_MAP] = keys_counts

    # Adjust
    df_norm_adj = df_norm.copy()
    
    _validate_merge(df_meta, df_norm_adj, mcol, ncol)
    return df_meta.merge(df_norm_adj, how='left', left_on=mcol, right_on=ncol)


def parse_g2pa_rexp(df_meta: pd.DataFrame, df_norm: pd.DataFrame) -> pd.DataFrame:
    df_meta = df_meta.copy()

    mcol = 'extractionID'
    ncol = 'sample_id'

    keys_counts = [f"g2-inc-pa-kallisto_counts-G2PA.{sid}.abundance.tsv" 
                   for sid in df_meta[mcol]]
    df_meta[COL_KALLISTO_MAP] = keys_counts

    # Adjust
    df_norm_adj = df_norm.copy()
    df_norm_adj.loc[df_norm_adj[ncol] == 'REXP1_LFeA',ncol] = 'REXP1LFeA'
    
    _validate_merge(df_meta, df_norm_adj, mcol, ncol)
    return df_meta.merge(df_norm_adj, how='left', left_on=mcol, right_on=ncol)


def parse_g2pa_st(df_meta: pd.DataFrame, df_norm: pd.DataFrame) -> pd.DataFrame:
    df_meta = df_meta.copy()

    prefix = 'g2-st-am-pa'
    regex = rf'(?<={prefix}\.)S(\d+)\.(C\d+)\.(.+um)\.(\w)'
    
    ignore = ['G2PA.S18C1.15m.0_2um.C']

    # Target regex pattern
    def formatter(match):
        s, c, f, r = match.groups()
        s = s.zfill(2)
        f_underscore = f.replace('.', '_')        
        k = f'G2PA.S{s}{c}.15m.{f_underscore}.{r}'
        kn = k
        return k, kn

    return _parse_station(
        df_meta, df_norm, regex, 
        col_norm='sample_name', key_formatter=formatter, 
        ignore_list=ignore)

def parse_g3ns_uw(df_meta: pd.DataFrame, df_norm: pd.DataFrame, extra_file: str) -> pd.DataFrame:
    df_meta = df_meta.copy()
    df_extra = pd.read_csv(extra_file)

    mcol = 'norm_key'
    ncol = 'sample_name'

    # Target regex pattern
    
    keys_counts = []
    keys_norm = []
    
    # Use zip to iterate over both the SampleID and your Depth column simultaneously
    # Replace 'Depth_m' with the exact name of your depth column
    for i, row in df_meta.iterrows():
        s, d, f, r = [row[col] for col in ['ID','Depth.m','Filter.um','Replicate']]
        s_unalias = df_extra.loc[df_extra['Alias2'] == s, 'Alias1'].values[0]
        s_splt = s_unalias.split()
        if len(s_splt) == 3:
            s1, s2, _ = s_splt
            s2 = s2.lstrip('#')
            s_sam = f'{s1}_{s2}'
        elif len(s_splt) == 2:
            s_sam, _ = s_splt
            s_sam += '_1'
        else:
            raise ValueError(f'Alias "{s_splt}" for sample "{s}" cannot be parsed')

        # Construct the count key dynamically
        f_und = str(f).rstrip('.0').replace('.', '_')
        kn = f'G3.UW.NS.{s_sam}.{d}m.{f_und}um.{r}'
        k = f'{kn}.tsv'
        keys_counts.append(k)
        keys_norm.append(kn)
        
    df_meta[COL_KALLISTO_MAP] = keys_counts
    df_meta[mcol] = keys_norm
    
    # Adjust
    df_norm_adj = df_norm.copy()
    # df_norm_adj.loc[df_norm_adj[ncol] == 'BD63','Sample_ID'] = 'BD60'
    
    _validate_merge(df_meta, df_norm_adj, mcol, ncol)
        
    return df_meta.merge(df_norm_adj, how='left', left_on=mcol, right_on=ncol)

def parse_g3pa_uw(df_meta: pd.DataFrame, df_norm: pd.DataFrame) -> pd.DataFrame:
    df_meta = df_meta.copy()

    mcol = 'norm_key'
    ncol = 'sample_name'

    # Target regex pattern
    
    keys_counts = []
    keys_norm = []
    
    for i, row in df_meta.iterrows():
        si, se, d, f, r = [row[col] for col in ['ID','extractionID','Depth.m','Filter.um','Replicate']]
    #     s_splt = si.split()
    #     if len(s_splt) == 4:
    #         s1, s2, _, _ = s_splt
    #         s2 = s2.lstrip('#')
    #         s_sam = f'{s1}_{s2}'
    #     elif len(s_splt) == 3:
    #         s_sam, _, _ = s_splt
    #         s_sam += '_1'
    #     else:
    #         raise ValueError(f'Alias "{s_splt}" for sample "{s}" cannot be parsed')

    #     # Construct the count key dynamically
    #     f_und = str(f).rstrip('.0').replace('.', '_')
    #     kn = f'G3.UW.PA.{s_sam}.{d}m.{f_und}um.{r}'

        k = f'G3PA.{se}.unstranded.abundance.tsv.gz'
        kn = f'G3PA.{se}.flash'
        keys_counts.append(k)
        keys_norm.append(kn)
        
    df_meta[COL_KALLISTO_MAP] = keys_counts
    df_meta[mcol] = keys_norm
    
    # Adjust
    df_norm_adj = df_norm.copy()
    # df_norm_adj.loc[df_norm_adj[ncol] == 'BD63','Sample_ID'] = 'BD60'
    
    _validate_merge(df_meta, df_norm_adj, mcol, ncol)
        
    return df_meta.merge(df_norm_adj, how='left', left_on=mcol, right_on=ncol)


def parse_g3papm_uw(df_meta: pd.DataFrame, df_norm: pd.DataFrame) -> pd.DataFrame:
    df_meta = df_meta.copy()

    mcol = 'norm_key'
    ncol = 'sample_name'

    # Target regex pattern
    keys_counts = []
    keys_norm = []
    
    for i, row in df_meta.iterrows():
        si, se, d, f, r = [row[col] for col in ['ID','extractionID','Depth.m','Filter.um','Replicate']]
        f = str(f).rstrip('.0')
        s_splt = si.split()
        if len(s_splt) == 3:
            s1, s2, _ = s_splt
            s2 = s2.lstrip('#')
            s_sam = f'{s1}_{s2}'
        elif len(s_splt) == 2:
            s_sam, _ = s_splt
            s_sam += '_1'
        else:
            raise ValueError(f'Alias "{s_splt}" for sample "{si}" cannot be parsed')

        base = f'G3.UW.PA.{s_sam}.{d}m_PM.{f}um.{r}'
        k = f'{base}.unstranded.abundance.tsv.gz'
        kn = f'{base}.flash'
        keys_counts.append(k)
        keys_norm.append(kn)
        
    df_meta[COL_KALLISTO_MAP] = keys_counts
    df_meta[mcol] = keys_norm
    
    # Adjust
    df_norm_adj = df_norm.copy()
    # df_norm_adj.loc[df_norm_adj[ncol] == 'BD63','Sample_ID'] = 'BD60'
    
    _validate_merge(df_meta, df_norm_adj, mcol, ncol)
        
    return df_meta.merge(df_norm_adj, how='left', left_on=mcol, right_on=ncol)


def parse_g3diel_uw(df_meta: pd.DataFrame, df_norm: pd.DataFrame) -> pd.DataFrame:
    df_meta = df_meta.copy()

    mcol = 'norm_key'
    ncol = 'sample_name'

    # Target regex pattern
    keys_counts = []
    keys_norm = []
    
    for i, row in df_meta.iterrows():
        si, se, d, f, r = [row[col] for col in ['ID','extractionID','Depth.m','Filter.um','Replicate']]
        s, _ = si.split()

        base = f'G3PA.diel.{s}.{r}'
        k = f'{base}.unstranded.abundance.tsv.gz'
        kn = f'{base}.flash'
        keys_counts.append(k)
        keys_norm.append(kn)
        
    df_meta[COL_KALLISTO_MAP] = keys_counts
    df_meta[mcol] = keys_norm
    
    # Adjust
    df_norm_adj = df_norm.copy()
    # df_norm_adj.loc[df_norm_adj[ncol] == 'BD63','Sample_ID'] = 'BD60'
    
    _validate_merge(df_meta, df_norm_adj, mcol, ncol)
        
    return df_meta.merge(df_norm_adj, how='left', left_on=mcol, right_on=ncol)


def parse_d1pa(df_meta: pd.DataFrame, df_norm: pd.DataFrame) -> pd.DataFrame:
    df_meta = df_meta.copy()

    mcol = 'norm_key'
    ncol = 'sample'

    # Target regex pattern
    keys_counts = []
    keys_norm = []
    
    for i, row in df_meta.iterrows():
        si, se, d, f, r = [row[col] for col in ['ID','extractionID','Depth.m','Filter.um','Replicate']]
        s, _ = si.split()

        k = f'D1PA.S{se}.abundance.tsv'
        kn = f'{se}'
        keys_counts.append(k)
        keys_norm.append(kn)
        
    df_meta[COL_KALLISTO_MAP] = keys_counts
    df_meta[mcol] = keys_norm
    
    # Adjust
    df_norm_adj = df_norm.copy()
    # df_norm_adj.loc[df_norm_adj[ncol] == 'BD63','Sample_ID'] = 'BD60'
    
    _validate_merge(df_meta, df_norm_adj, mcol, ncol)
        
    return df_meta.merge(df_norm_adj, how='left', left_on=mcol, right_on=ncol)

# -----------------------------------------------------------------------------
# Function Registry
# -----------------------------------------------------------------------------
# This maps your string keys to the actual Python functions
FUNCTION_REGISTRY = {
    "g1-st-am-ns": parse_g1ns,
    "g1-st-am-pa": parse_g1pa,
    "g2-st-am-ns": parse_g2ns_st,
    "g2-st-am-pa": parse_g2pa_st,
    "g2-ctd-ns": parse_g2ns_dcm,
    "g2-ctd-pa": parse_g2pa_dcm,
    "g2-inc-ns": parse_g2ns_rexp,
    "g2-inc-pa": parse_g2pa_rexp,
    "g3-uw-am-ns": parse_g3ns_uw,
    "g3-uw-am-pa": parse_g3pa_uw,
    "g3-uw-pm-pa": parse_g3papm_uw,
    "g3-diel-pa": parse_g3diel_uw,
    "d1-st-pa": parse_d1pa,
    # "G1NS": parse_g1ns,
    # "G1PA": parse_g1pa,
    # "G2NS.ST": parse_g2ns_st,
    # "G2PA.ST": parse_g2pa_st,
    # "G2NS.DCM": parse_g2ns_dcm,
    # "G2PA.DCM": parse_g2pa_dcm,
    # "G2NS.REXP": parse_g2ns_rexp,
    # "G2PA.REXP": parse_g2pa_rexp,
    # "G3NS.UW": parse_g3ns_uw,
    # "G3PA.UW": parse_g3pa_uw,
}
COL_KALLISTO_MAP = 'kallisto_fn'
# -----------------------------------------------------------------------------
# Click CLI Command Configuration
# -----------------------------------------------------------------------------
@click.command()
@click.option(
    "--meta-file",
    "-m",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
    help="Path to the CSV file containing metadata.",
)
@click.option(
    "--norm-file",
    "-n",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
    help="Path to the CSV file containing normalization factors.",
)
@click.option(
    "--counts-col-names-file",
    "-c",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
    help="Path to a file with the list of column names to match the metadata with.",
)
@click.option(
    "--extra-file",
    "-e",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="Path to a file with extra info to pass in parsing filenames.",
)

@click.option(
    "--method-key",
    "-k",
    type=click.Choice(list(FUNCTION_REGISTRY.keys()), case_sensitive=True),
    required=True,
    help="Key deciding which parsing function to execute.",
)
@click.option(
    "--output-file",
    "-o",
    type=click.Path(dir_okay=False, writable=True),
    required=True,
    help="Path where the output CSV will be saved.",
)

def main(meta_file, norm_file, counts_col_names_file, extra_file, method_key, output_file):
    """
    Load metadata and normalization CSVs, apply a selected parsing function
    based on a method key, and save the resulting DataFrame to an output file.
    """
    # try:
    # 1. Load the CSVs into pandas DataFrames
    click.echo(f"Loading files:\n - Metadata: {meta_file}\n - Norm Factors: {norm_file}")
    df_meta = pd.read_csv(meta_file)
    df_norm = pd.read_csv(norm_file)

    # 2. Select the correct function using the string key
    # (click.Choice already validates that the key exists in our dict)
    parse_function = FUNCTION_REGISTRY[method_key]

    # 3. Execute the function
    if extra_file is None:
        result_df = parse_function(df_meta, df_norm)
    else:
        result_df = parse_function(df_meta, df_norm, extra_file=extra_file)
    result_df['batch_meta'] = method_key
    ## Check
    # Names from the metadata tables to ignore when mapping to kallisto filenames
    ignore = [
        'REXP1LFeA','G2PA.S18C1.15m.0_2um.C','G3PA.UW6.unstranded.abundance.tsv.gz', 
        'G3PA.UW58.unstranded.abundance.tsv.gz', 'G3PA.UW18.unstranded.abundance.tsv.gz',
        'G3PA.diel.S4C4.C.unstranded.abundance.tsv.gz', 
        'G3PA.diel.S4C13.A.unstranded.abundance.tsv.gz', 
        'G3PA.diel.S4C21.A.unstranded.abundance.tsv.gz', 
        'G3PA.diel.S4C8.A.unstranded.abundance.tsv.gz',
    ]

    df_count_cols = pd.read_csv(counts_col_names_file, header=None)
    df_count_cols['colname'] = df_count_cols[0]
    _validate_merge(result_df, df_count_cols, COL_KALLISTO_MAP, 'colname',ignore_list=ignore)


    # 4. Write the new dataframe to the output filename
    click.echo(f"Saving results to: {output_file}")
    result_df.to_csv(output_file, index=False)
    
    click.echo("Processing complete successfully!")

    # except Exception as e:
    #     click.echo(f"Error occurred during execution: {e}", err=True)
    #     sys.exit(1)


if __name__ == "__main__":
    main()