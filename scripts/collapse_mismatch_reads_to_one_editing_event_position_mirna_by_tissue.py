import pandas as pd
import sys
import os


"""
Load files from command line arguments
"""


input_file = sys.argv[1]
out_file = sys.argv[2]

df = pd.read_csv(input_file, sep='\t')


"""
Define functions
"""


def split_mismatch_name(mismatch):
    if str(mismatch) == 'nan':
        mismatch = 'None_None'
    mismatch = mismatch.split('_')
    return mismatch[0], str(mismatch[1])


def split_edit_from_to(mismatch):
    edit_from = mismatch.split('to')[0]
    edit_to = mismatch.split('to')[1].split('_')[0]
    return edit_from, edit_to


"""
Collapse the mismatch reads files into one editing 
event / position per miRNA
"""


# Sum total number of reads per miRNA
df['total'] = df.groupby(['miRNA','tissue'])['nr_reads'].transform('sum')

# Transform dataframe to long format
df = df.melt(id_vars=['ID', 'nr_reads', 'miRNA', 'arm', 'mature_star', 
                     'species', 'tissue', 'read_start', 'read_seq', 'total'],
                     value_vars=['from_to_1', 'from_to_2', 'from_to_3'],
                     value_name='mismatch').drop('variable', axis=1)

# Sum total number for each mismatch per miRNA position
total_mismatch = df.groupby(['miRNA', 'mismatch'])['nr_reads'].transform('sum').rename('total_mismatch').reset_index()
df['total_mismatch'] = total_mismatch['total_mismatch']
df = df[['miRNA', 'arm', 'mature_star', 'species', 'tissue', 'total', 'mismatch', 'total_mismatch']].drop_duplicates().dropna().sort_values(by=['total', 'total_mismatch'], ascending=False)

# Edit the columns for easier downstream analysis
df = df.reset_index().drop('index', axis=1)
df['mismatch_pos'] = pd.Series(map(lambda x: split_mismatch_name(x)[1], df['mismatch']))
df['edit_from'] = pd.Series(map(lambda x: split_edit_from_to(x)[0], df['mismatch']))
df['edit_to'] = pd.Series(map(lambda x: split_edit_from_to(x)[1], df['mismatch']))
df['mismatch'] = pd.Series(map(lambda x: split_mismatch_name(x)[0], df['mismatch']))
df = df[['miRNA', 'arm', 'mature_star', 'species', 'tissue', 'total', 'mismatch', 'edit_from', 'edit_to', 'mismatch_pos', 'total_mismatch']]

# Write .tsv file output
df.to_csv(out_file, sep='\t')