import argparse
import json
import re
from Bio import SeqIO
import pandas as pd
import numpy as np


def parse_arguments():
    parser = argparse.ArgumentParser(description="Description of your script")
    parser.add_argument("species_dict", help="Path to the species dictionary file")
    parser.add_argument("ref", help="Path to the reference fasta file")
    parser.add_argument("sam_file", help="Path to the sam file")
    parser.add_argument("out_file", help="Path to the output csv file")
    parser.add_argument("species_id", help="Species ID")
    parser.add_argument("tissue_name", help="Tissue name")
    return parser.parse_args()


def load_species_dict(species_dict_path):
    with open(species_dict_path, 'r') as f:
        return json.load(f)


def fetch_mismatch(read_start, y, mismatch_pos, n, mirna, sequence, ref_seq):
    canon_seq = ref_seq[mirna].seq[5:-5]

    k = np.array([[2, 4, 6], [1, 3, 5]])

    canon_pos = np.zeros(3)

    for i in range(len(y)):
        mis_nuc = mismatch_pos[k[0, i]]
        mis_pos = 1 + int(mismatch_pos[k[1, i]]) + (int(y[i-1][1]) if i > 0 else read_start - 1)
        if mis_pos > len(canon_seq):
            mis_pos = '+' + str(mis_pos - len(canon_seq))
        y[i] = [mis_nuc, str(mis_pos)]
        if i > 0:
            canon_pos[i] = 1 + int(mismatch_pos[k[1, i]]) + canon_pos[i-1]

    if len(y) < n:
        return ''
    else:
        last_mismatch = y[n-1]
        last_mis_pos = int(last_mismatch[1])
        adjusted_pos = last_mis_pos - read_start + 1
        return f"{last_mismatch[0]}to{sequence[adjusted_pos]}_{last_mismatch[1]}"


def parse_sam(file_path, species, tissue):
    # Implementation of parse_sam function
    skip_lines = 0
    with open(file_path, 'r') as sam_file:
        for line in sam_file:
            if line.startswith('@HD') or line.startswith('@SQ') or line.startswith('@PG'):
                skip_lines += 1
            else:
                break
    sam_df = pd.read_csv(file_path, sep='\t', skiprows=skip_lines, header=None)
    sam_df = sam_df.rename(columns={0: 'ID', 2: 'miRNA', 3: 'read_start', 9: 'read_seq', 12: 'MD:mis_pos', 13: 'NM:i:#'})            
    sam_df['read_start'] = sam_df['read_start'] - 5  # Adjust for 5prime +5nt extension reference
    sam_df['species'] = species
    sam_df['tissue'] = tissue
    sam_df['nr_reads'] = pd.Series(map(lambda x: split_ID(x), sam_df['ID']))
    sam_df['arm'] = pd.Series(map(lambda x: x.split('_')[1], sam_df['miRNA']))
    sam_df['mature_star'] = pd.Series(map(lambda x: mature_star(x), sam_df['arm']))
    sam_df['nr_mismatches'] = pd.Series(map(lambda x: int(x.split(':')[2]), sam_df['NM:i:#']))
    sam_df['mismatch_pos'] = pd.Series(map(lambda x: re.split('(\d+)', x), sam_df['MD:mis_pos']))
    sam_df['y'] = pd.Series(map(lambda x: [0] * x, sam_df['nr_mismatches']))
    sam_df['from_to_1'] = pd.Series(map(lambda read_start, y, mismatch_pos, mirna, sequence: fetch_mismatch(read_start, y, mismatch_pos, 1, mirna, sequence), sam_df['read_start'], sam_df['y'], sam_df['mismatch_pos'], sam_df['miRNA'], sam_df['read_seq']))
    sam_df['from_to_2'] = pd.Series(map(lambda read_start, y, mismatch_pos, mirna, sequence: fetch_mismatch(read_start, y, mismatch_pos, 2, mirna, sequence), sam_df['read_start'], sam_df['y'], sam_df['mismatch_pos'], sam_df['miRNA'], sam_df['read_seq']))
    sam_df['from_to_3'] = pd.Series(map(lambda read_start, y, mismatch_pos, mirna, sequence: fetch_mismatch(read_start, y, mismatch_pos, 3, mirna, sequence), sam_df['read_start'], sam_df['y'], sam_df['mismatch_pos'], sam_df['miRNA'], sam_df['read_seq']))
    sam_df = sam_df[['ID', 'nr_reads', 'miRNA', 'arm', 'mature_star', 'species', 'tissue', 'read_start', 'read_seq', 'from_to_1', 'from_to_2', 'from_to_3']]
    return sam_df


def count_total_reads(mirna, counts):
    # Implementation of count_total_reads function
    total = pd.Series([0]*len(mirna))
    local_df = pd.DataFrame([mirna, counts, total], index=['miRNA', 'counts', 'total']).T

    for x in local_df.miRNA.unique():
        local_df.loc[local_df.miRNA == x, 'total'] = local_df[local_df.miRNA == x].counts.sum()
    return local_df.total


def main():
    args = parse_arguments()
    species_names = load_species_dict(args.species_dict)
    ref_seq = SeqIO.to_dict(SeqIO.parse(args.ref, 'fasta'))
    df = parse_sam(args.sam_file, species_names[args.species_id.capitalize()], args.tissue_name)
    df.to_csv(args.out_file, sep='\t')

if __name__ == "__main__":
    main()
