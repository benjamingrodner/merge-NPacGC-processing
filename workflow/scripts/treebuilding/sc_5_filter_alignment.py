#!/usr/bin/env python3

import argparse

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


def main():
    # parser = argparse.ArgumentParser()
    # parser.add_argument("fn_aln", type=str, help="alignment file")
    # parser.add_argument("thresh_len", type=int, help="minimum length")
    # args = parser.parse_args()

    fn_aln = 'mafft/fbp1_btuF_desE_hmuT_fecB_btuF2_mtsA_clustered_source_tax_subdb1pct-1esz_2chu_3eiw_2r79_1n2z/fbp1_btuF_desE_hmuT_fecB_btuF2_mtsA_clustered_source_tax_subdb1pct-1esz_2chu_3eiw_2r79_1n2z.faa.aln.known_trim.clipkit'
    thresh_len = 150
    fn_trim = f'{fn_aln}.minlen{thresh_len}'
    fn_tooshort = f'{fn_aln}.minlen{thresh_len}_tooshort'
    dict_header_seq = fasta_to_dict(fn_aln)
    dict_out = {}
    tooshorts = []
    for h, s in dict_header_seq.items():
        ngaps = s.count('-')
        slen = len(s) - ngaps
        if slen > thresh_len:
            dict_out[h] = s
        else:
            tooshorts.append(h)
    with open(fn_trim, 'w') as f:
        for h, s in dict_out.items():
            f.write(f'>{h}\n{s}\n')
    with open(fn_tooshort, 'w') as f:
        for h in tooshorts:
            f.write(f'{h}\n')
    return

if __name__ == "__main__":
    main()