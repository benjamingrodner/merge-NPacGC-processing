
## Prompt to gemini:
# You are an expert Python data engineer with strong experience in ibis, 
# DuckDB, parquet, and parsing HMMER output. Write production-quality Python code to solve the following task.

# Task
# I have:

# a large parquet table with columns including target_id, KO, and optionally E-value and/or score;
# several HMMER hmmsearch --tblout result tables, each associated with a specific gene name.
# For each row in the parquet table:

# if target_id appears in any HMMER result table, consider the corresponding gene name as a candidate replacement for KO;
# replace KO with that gene name only if the HMMER hit is better than the existing row according to the relevant metric:
# lower E-value is better;
# higher score is better.
# The code must use ibis with a DuckDB backend to read and transform the parquet table.

# Requirements
# Read the large parquet table lazily with ibis + DuckDB.
# Parse multiple HMMER --tblout files, each associated with a gene name.
# Normalize the HMMER tables into a consistent schema.
# Combine all HMMER hits into one table.
# If multiple HMMER rows match the same target_id, reduce them to the best hit per target using the correct comparison rule.
# Join the HMMER hits to the parquet table on target_id.
# Update KO only when the HMMER hit is better than the existing row.
# Preserve all original rows and columns from the parquet table.
# Write the updated result to a parquet file.
# Include clear, minimal comments and make the code runnable as a script.
# Important details
# Assume the HMMER --tblout files are tabular but may need parsing because they are whitespace-delimited rather than clean TSV.
# Handle missing values robustly.
# The code should be efficient for large tables and avoid loading the full parquet table into pandas unless absolutely necessary.
# Prefer clean ibis expressions over manual row-by-row Python loops.
# If E-value and score both exist, implement a sensible rule and make it explicit in the code.
# If the update rule is ambiguous, choose a clear precedence and document it briefly in a comment.
# Include command-line arguments for:
# input parquet path,
# output parquet path,
# one or more HMMER --tblout files,
# gene name for each file, or a way to infer it from filename.
# Suggested output structure
# Return:

# A short explanation of the approach.
# A complete Python script.
# Any assumptions made.
# Code quality requirements
# Use functions for parsing and transformation.
# Include type hints where helpful.
# Validate inputs and give informative errors.
# Keep the implementation practical and concise.
# Do not use placeholder pseudocode; produce real code.
# Clarifying assumptions to make if needed
# If the input parquet has both E-value and score, the code should define whether:

# either metric can trigger replacement,
# both must agree, or
# one metric takes precedence.
# Choose a reasonable approach and implement it consistently.
# Use ibis for the main join/update logic, and only use pandas for parsing the small HMMER tables if necessary.

### Code edited by perplexity

### Edits by Ben grodner

import argparse
import sys
from pathlib import Path
from typing import List, Optional

import ibis
import pandas as pd

JOIN_KEY = config['keys_kallisto']['join']
COL_6TR = config['keys_kofam']['col_6tr']

def parse_hmmer_tblout(filepath: Path, gene_name: str) -> pd.DataFrame:
    if not filepath.exists():
        raise FileNotFoundError(f"HMMER file not found: {filepath}")

    df = pd.read_csv(
        filepath,
        sep=r"\s+",
        comment="#",
        header=None,
        engine="python"
    )
    columns = [
        "target_id_6tr_hmm", "target_accession", "query_name", "query_accession",
        "hmmer_evalue", "hmmer_score", "full_bias", "best_domain_evalue",
        "best_domain_score", "best_domain_bias", "exp", "reg", "clu", "ov",
        "env", "dom", "rep", "inc", "description"
    ]
    if not df.shape[1] == len(columns):
        raise ValueError(f"Hmmer table in {filepath} does not have the correct number of columns. Needed: {len(columns)}, Present: {df.shape[1]}")
    
    df.columns = columns
    df['hmmer_gene'] = gene_name

    return df[["target_id_6tr_hmm", "hmmer_gene", "hmmer_evalue", "hmmer_score"]]


def parse_hmmer_inputs(
    tblout_files: List[str],
    gene_names: Optional[List[str]] = None,
) -> pd.DataFrame:
    dfs = []

    for idx, file_str in enumerate(tblout_files):
        file_path = Path(file_str)
        gene = gene_names[idx] if gene_names and idx < len(gene_names) else file_path.stem.split(".")[0]

        df = parse_hmmer_tblout(file_path, gene)
        if not df.empty:
            df[JOIN_KEY] = df["target_id_6tr_hmm"].str.replace(r'_\d+$','',regex=True)
            dfs.append(df)

    if not dfs:
        return pd.DataFrame(
            columns=[JOIN_KEY, "target_id_6tr_hmm", "hmmer_gene", "hmmer_evalue", "hmmer_score"]
        )

    return pd.concat(dfs, ignore_index=True)


def get_metric_columns(columns: List[str]):
    evalue_col = "eval_ko" if "eval_ko" in columns else None
    score_col = "score" if "score" in columns else None
    if (evalue_col is None) and (score_col is None):
        raise ValueError(f"Either 'eval_ko' or 'score' must be a column in the data table, neither is present.")
    return evalue_col, score_col


# def build_best_hmmer_expr(hmmer_raw):
#     return (
#         hmmer_raw
#         .mutate(
#             rn=ibis.row_number().over(
#                 ibis.window(
#                     group_by=[hmmer_raw[config['colname_6tr']]],
#                     order_by=[hmmer_raw.hmmer_evalue.asc(), hmmer_raw.hmmer_score.desc()],
#                 )
#             )
#         )
#         .filter(lambda t: t.rn == 1)
#         .drop("rn")
#     )

def get_best_hmmer_df(hmmer_raw):
    return (
        hmmer_raw
        .sort_values(
            by=[JOIN_KEY, 'hmmer_evalue', 'hmmer_score'],
            ascending=[True, True, False]
        )
        .groupby(JOIN_KEY, as_index=False)
        .head(1)
    )


def build_replacement_expr(joined, columns: List[str]):
    evalue_col, score_col = get_metric_columns(columns)

    hmmer_present = joined.hmmer_gene.notnull()
    
    new_score = None
    new_eval = None
    if evalue_col and score_col:
        orig_e = joined[evalue_col]
        orig_s = joined[score_col]
        better = (
            orig_e.isnull()
            | (joined.hmmer_score > orig_s)
            | ((joined.hmmer_score == orig_s) & (joined.hmmer_evalue > orig_e))
        )
        new_score = ibis.ifelse(
            hmmer_present & better, joined.hmmer_score, joined[score_col]
        )
        new_eval = ibis.ifelse(
            hmmer_present & better, joined.hmmer_evalue, joined[evalue_col]
        )
    elif evalue_col:
        orig_e = joined[evalue_col]
        better = orig_e.isnull() | (joined.hmmer_evalue < orig_e)
        new_eval = ibis.ifelse(
            hmmer_present & better, joined.hmmer_evalue, joined[evalue_col]
        )
    elif score_col:
        orig_s = joined[score_col]
        better = orig_s.isnull() | (joined.hmmer_score > orig_s)
        new_score = ibis.ifelse(
            hmmer_present & better, joined.hmmer_score, joined[score_col]
        )

    new_ko = ibis.ifelse(hmmer_present & better, joined.hmmer_gene, joined.KO)
    new_6tr = ibis.ifelse(hmmer_present & better, joined["target_id_6tr_hmm"], joined[COL_6TR])
    return new_ko, new_score, new_eval, new_6tr


def project_updated_table(joined, columns: List[str], new_ko, new_score, new_eval, new_6tr):
    projected = []
    for col in columns:
        if col == "KO":
            projected.append(new_ko.name("KO"))
        elif col == 'score':
            if new_score is None:
                raise ValueError("'score' colum present in main table, but replacement expression not created")
            projected.append(new_score.name("eval_score"))
        elif col == 'eval_ko':
            if new_eval is None:
                raise ValueError("'eval_ko' colum present in main table, but replacement expression not created")
            projected.append(new_eval.name("eval_ko"))
        elif col == COL_6TR:
            projected.append(new_6tr.name(COL_6TR))
        else:
            projected.append(joined[col])
    return joined.select(*projected)


def build_updated_table(con: ibis.BaseBackend, parquet_path: str, hmmer_df: pd.DataFrame):
    main_tbl = con.read_parquet(parquet_path)

    if hmmer_df.empty:
        print("No hmmer hits, no rows edited.")
        return main_tbl


    if "KO" not in main_tbl.columns:
        raise ValueError("Input parquet must contain a KO column.")
    if JOIN_KEY not in main_tbl.columns:
        raise ValueError(f"Input parquet must contain a {JOIN_KEY} column.")

    # hmmer_raw = con.create_table("hmmer_raw", hmmer_df, temp=True)
    best_hmmer = get_best_hmmer_df(hmmer_df)
    best_hmmer = con.create_table("hmmer_raw", best_hmmer, temp=True)


    joined = main_tbl.left_join(best_hmmer, JOIN_KEY)
    new_ko, new_score, new_eval, new_6tr = build_replacement_expr(
        joined, list(main_tbl.columns)
    )
    return project_updated_table(
        joined, list(main_tbl.columns), new_ko, new_score, new_eval, new_6tr
    )


rule annotate_custom_hmms:
    input:
        fn_merged_all = fmt_merged_all,
    output:
        fn_custom_ann = fmt_custom_ann,
    threads: 20
    resources:
        mem="100G"
    run:
        fns_tbl = list(config['dict_gene_fnhmm'].values())
        fns_gene = list(config['dict_gene_fnhmm'].keys())
        hmmer_df = parse_hmmer_inputs(fns_tbl, fns_gene)

        con = ibis.duckdb.connect(memory_limit=resources.mem, threads=threads)
        result_tbl = build_updated_table(con, input.fn_merged_all, hmmer_df)

        result_tbl.to_parquet(output.fn_custom_ann)



