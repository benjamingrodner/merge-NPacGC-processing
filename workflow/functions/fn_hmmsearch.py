import argparse
import json
import gzip
import re

def extract_sequences_from_fasta(fasta_path, header_names):
    """
    Extracts sequences from a gzipped or plain FASTA file by header names.

    Args:
        fasta_path (str): Path to the FASTA file (can be gzipped or plain).
        header_names (set or list): Sequence headers to extract (without '>' and possibly without whitespace).

    Returns:
        dict: Mapping from header name to sequence string.
    """
    def file_opener(path):
        if path.endswith('.gz'):
            return gzip.open(path, 'rt')
        else:
            return open(path, 'r')

    results = {}
    current_header = None
    current_seq = []
    keep = False
    wanted = set([re.sub(r'_\d$','',h) for h in header_names])  # remove 6tr id
    with file_opener(fasta_path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith('>'):
                if current_header and keep:
                    results[current_header] = ''.join(current_seq)
                header = line[1:].split()[0]
                header = re.sub(r'_\d$','',header)  # remove 6tr id
                keep = header in wanted
                current_header = header
                current_seq = []
            elif keep:
                current_seq.append(line)
        # Capture last entry if needed
        if current_header and keep:
            results[current_header] = ''.join(current_seq)
    return results


def parse_hmmsearch_tblout_summary(tblout_path):
    """
    Parses a HMMER hmmsearch/hmmscan --tblout summary file.

    Args:
        tblout_path (str): Path to the --tblout file.

    Returns:
        List[dict]: List of hits, each as a dictionary with keys:
            'target_name', 'target_accession', 
            'query_name', 'query_accession', 
            'seq_evalue', 'seq_score', 'seq_bias', 
            'dom_evalue', 'dom_score', 'dom_bias', 
            'exp', 'reg', 'clu', 'ov', 'env', 'dom', 
            'rep', 'inc', 'description'
    """
    hits = []
    with open(tblout_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue  # skip header/comments
            parts = line.strip().split()
            if len(parts) < 18:
                print(f"Warning: line with unexpected format ({len(parts)} cols):")
                print(parts)
                continue

            hit = {
                'target_name': parts[0],
                'target_accession': parts[1],
                'query_name': parts[2],
                'query_accession': parts[3],
                'seq_evalue': float(parts[4]),
                'seq_score': float(parts[5]),
                'seq_bias': float(parts[6]),
                'dom_evalue': float(parts[7]),
                'dom_score': float(parts[8]),
                'dom_bias': float(parts[9]),
                'exp': float(parts[10]),
                'reg': int(parts[11]),
                'clu': int(parts[12]),
                'ov': int(parts[13]),
                'env': int(parts[14]),
                'dom': int(parts[15]),
                'rep': int(parts[16]),
                'inc': int(parts[17]),
                'description': ' '.join(parts[18:]) if len(parts) > 18 else ''
            }
            hits.append(hit)
    return hits


def extract_fasta_headers(fasta_file):
    """
    Extracts headers from a FASTA file.

    Parameters:
        fasta_file (str): Path to the FASTA file.

    Returns:
        list: A list of headers without the leading '>'.
    """
    headers = []
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                headers.append(line[1:].strip())  # Remove '>' and any extra spaces
    return headers


def fasta_to_dict(fasta_file):
    """
    Reads a FASTA file and returns a dictionary of {header: sequence}.

    Parameters:
        fasta_file (str): Path to the FASTA file.

    Returns:
        dict: Keys are headers (without '>'), values are sequences as strings.
    """
    fasta_dict = {}
    header = None
    seq_chunks = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue  # skip empty lines
            if line.startswith(">"):
                # Save the previous entry
                if header:
                    fasta_dict[header] = ''.join(seq_chunks)
                header = line[1:].strip()  # remove '>'
                seq_chunks = []
            else:
                seq_chunks.append(line)

        # Don't forget the last entry
        if header:
            fasta_dict[header] = ''.join(seq_chunks)

    return fasta_dict