import sys
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import re

# Make dictionary of species names and abreviations used in MirGeneDB
species_names = {
 "Hsa" : "Human", 
 "Mmu" : "Mouse", 
 "Rno" : "Rat", 
 "Cpo" : "Guinea_Pig", 
 "Ocu" : "Rabbit", 
 "Dno" : "Armadillo", 
 "Gga" : "Chicken", 
 "Dre" : "Zebrafish", 
 "Mml" : "Rhesus_monkey", 
 "Cfa" : "Dog", 
 "Bta" : "Cow", 
 "Ete" : "Lesser_hedgehog_tenrec", 
 "Cli" : "Rock_pigeon", 
 "Ami" : "American_aligator", 
 "Cpi" : "Western_painted_turtle", 
 "Aca" : "Green_anole_lizard", 
 "Xtr" : "Tropical_clawed_frog", 
 "Bfl" : "Florida_lancet", 
 "Sko" : "Acorn_worm", 
 "Spu" : "Purple_sea_urchin", 
 "Dme" : "Fruit_fly_Dme", 
 "Dmo" : "Fruit_flu_Dmo",
 "Dan" : "Fruit_fly_Dan",
 "Dpu" : "Common_water_flea", 
 "Isc" : "Black-legged_tick", 
 "Cel" : "Roundworm_Cel", 
 "Lgi" : "Owl_limpet", 
 "Cgi" : "Pacific_oyster", 
 "Cte" : "Polychaete_worm", 
 "Efe" : "Common_brandling_worm", 
 "Pmi" : "Bar_starfish", 
 "Tca" : "Red flour beetle", 
 "Asu" : "Large roundworm", 
 "Pfl" : "Acorn worm 2", 
 "Mdo" : "Opossum", 
 "Sha" : "Tasmanian devil",
 "Aae" : "Yellow fever mosquito",
 "Hme" : "Longwing",
 "Bge" : "Cockroach",
 "Cbr" : "Roundworm_Cbr",
 "Lan" : "Lingula",
 "Cin" : "Sea Squirt",
 "Tgu" : "Zebra finch",
 "Oan" : "Platypus",
 "Sto" : "Cloudy Catshark",
 "Aqu" : "Amphimedon",
 "Bla" : "European lancet",
 "Cmi" : "Australian ghostshark",
 "Dsi" : "Fruit_fly_Dsi",
 "Dya" : "Fruit_fly_Dya",
 "Ebu" : "Inshore hagfish",
 "Emu" : "Muellers freshwater sponge",
 "Gja" : "Schlegels Japanese gecko",
 "Gmo" : "Cod",
 "Hvu" : "Freshwater-polyp",
 "Lch" : "Coelacanth",
 "Loc" : "Spotted gar",
 "Lpo" : "Atlantic horseshoe crab",
 "Mal" : "Asian swamp eel",
 "Mun" : "Microcaecelia",
 "Nve" : "Starlet sea anemone",
 "Pbv" : "Burmese python",
 "Pma" : "Sea Lamprey",
 "Spt" : "Tuatara",
 "Tni" : "Pufferfish",
 "Xla" : "African clawed Frog",
 "Csc" : "Arizona bark skorpion",
 "Obi" : "California two-spot octopus",
 "Ovu" : "Common octopus",
 "Bpl" : "Rotifer",
 "Dma" : "Large commom water flea",
 "Esc" : "Hawaiian bobtail squid",
 "Npo" : "Chambered Nautilus",
 "Sme" : "Freshwater planarian",
 "Xbo" : "Xenoturbella"
}

reference_seq_dict = {}

for file in os.listdir('mirgenedb_merged_extended_genome/'):
    if file.endswith('.fas'):
        record_dict  = SeqIO.to_dict(SeqIO.parse('mirgenedb_merged_extended_genome/' + file, 'fasta'))
        species_name = species_names[file.split('.')[0].capitalize()]
        reference_seq_dict[species_name] = record_dict 

def fetch_mismatch(read_start, y, mismatch_pos, n, mirna, species, sequence):
    canon_seq = reference_seq_dict[species][mirna].seq[5:-5]
    read_seq  = sequence
    if len(read_seq) > len(canon_seq):
        overhang_len = len(read_seq) - len(canon_seq)
    else:
        overhang_len = 0

    k = np.matrix([[2, 4, 6], [1, 3, 5]])

    canon_pos = [0]*3

    for i in range(len(y)):
        if i > 0:
            mis_nuc = mismatch_pos[k[0, i]]
            mis_pos = 1 + int(mismatch_pos[k[1, i]]) + int(y[i-1][1])
            if mis_pos > len(canon_seq):
                mis_pos = '+' + str(mis_pos - len(canon_seq))
            y[i] = [ mis_nuc, str(mis_pos)]
            if i == 1:
                canon_pos[i] = 1 + int(mismatch_pos[k[1,i]]) + canon_pos[i-1]
            if i == 2:
                canon_pos[i] = 1 + int(mismatch_pos[k[1,i]]) + canon_pos[i-1]
        else:
            mis_nuc = mismatch_pos[k[0, i]]
            mis_pos = 1 + int(mismatch_pos[k[1, i]]) + read_start - 1
            if mis_pos > len(canon_seq):
                mis_pos = '+' + str(mis_pos - len(canon_seq))
            y[0] = [ mis_nuc, str(mis_pos)]
            canon_pos[i] = int(mismatch_pos[k[1, i]]) + read_start - 1

    if len(y) < n:
        return ''
    else:
        return str(y[n-1][0]) + 'to' + str(sequence[canon_pos[n-1] - read_start+1]) + '_' + str(y[n-1][1])

def split_ID(ID):
    ID = ID.replace('#', '-').replace('x', '-').split('-')[1]
    return int(ID)

def mature_star(arm):
    if arm.endswith('p'):
        return 'mature'
    if arm.endswith('*'):
        return 'star'
    else:
        return 'NA'

def parse_sam(file, species, tissue):
    skip_lines = 0
    with open(file, 'r') as sam_file:
        for line in sam_file:
            if line.startswith('@HD') or line.startswith('@SQ') or line.startswith('@PG'):
                skip_lines += 1
            else:
                break
    sam_df = pd.read_csv(file, sep = '\t', skiprows=skip_lines, header=None)
    sam_df = sam_df.rename(columns={0: 'ID', 2: 'miRNA', 3: 'read_start', 9: 'read_seq', 12: 'MD:mis_pos', 13: 'NM:i:#'})            
    sam_df['read_start']    = sam_df['read_start'] - 5 # Adjust for 5prime +5nt extension reference
    sam_df['species']       = species
    sam_df['tissue']        = tissue
    sam_df['nr_reads']      = pd.Series(map(lambda x : split_ID(x), sam_df['ID']))
    sam_df['arm']           = pd.Series(map(lambda x : x.split('_')[1], sam_df['miRNA']))
    sam_df['mature_star']   = pd.Series(map(lambda x : mature_star(x), sam_df['arm'] ))
    sam_df['nr_mismatches'] = pd.Series(map(lambda x : int(x.split(':')[2]), sam_df['NM:i:#']))
    sam_df['mismatch_pos']  = pd.Series(map(lambda x : re.split('(\d+)', x), sam_df['MD:mis_pos']))
    sam_df['y']             = pd.Series(map(lambda x : [0] * x, sam_df['nr_mismatches']))
    sam_df['from_to_1']     = pd.Series(map(lambda read_start, y, mismatch_pos, mirna, sequence: fetch_mismatch(read_start, y, mismatch_pos, 1, mirna, species, sequence), sam_df['read_start'], sam_df['y'], sam_df['mismatch_pos'], sam_df['miRNA'], sam_df['read_seq']))
    sam_df['from_to_2']     = pd.Series(map(lambda read_start, y, mismatch_pos, mirna, sequence: fetch_mismatch(read_start, y, mismatch_pos, 2, mirna, species, sequence), sam_df['read_start'], sam_df['y'], sam_df['mismatch_pos'], sam_df['miRNA'], sam_df['read_seq']))
    sam_df['from_to_3']     = pd.Series(map(lambda read_start, y, mismatch_pos, mirna, sequence: fetch_mismatch(read_start, y, mismatch_pos, 3, mirna, species, sequence), sam_df['read_start'], sam_df['y'], sam_df['mismatch_pos'], sam_df['miRNA'], sam_df['read_seq']))
    sam_df = sam_df[['ID', 'nr_reads', 'miRNA', 'arm', 'mature_star', 'species', 'tissue', 'read_start', 'read_seq', 'from_to_1', 'from_to_2', 'from_to_3']]
    return(sam_df)

def count_total_reads(mirna, counts):
    total = pd.Series([0]*len(mirna))
    local_df = pd.DataFrame([mirna, counts, total], index=['miRNA', 'counts', 'total']).T

    for x in local_df.miRNA.unique():
        local_df.loc[local_df.miRNA == x, 'total'] = local_df[local_df.miRNA==x].counts.sum()
    return local_df.total

def split_mismatch_name(mismatch):
    if str(mismatch) == 'nan':
        mismatch = 'None_None'
    mismatch = mismatch.split('_')
    return mismatch[0], str(mismatch[1])

def split_edit_from_to(mismatch):
    edit_from = mismatch.split('to')[0]
    edit_to   = mismatch.split('to')[1].split('_')[0]
    return edit_from, edit_to


# Loop through species directories, and merge all samfiles into one pandas dataframe, written to .csv
print('Counting Mismatches')

for directory in os.listdir('mirgenedb_merged_extended_sam_filtered_no_outside_pos2_18_mismatch/'):
    if not os.path.isdir('mirgenedb_merged_extended_sam_filtered_no_outside_pos2_18_mismatch/' + str(directory)):
        continue
    print('nr_samfiles: ', len(os.listdir('mirgenedb_merged_extended_sam_filtered_no_outside_pos2_18_mismatch/' + str(directory))))
    if len(os.listdir('mirgenedb_merged_extended_sam_filtered_no_outside_pos2_18_mismatch/' + str(directory))) == 0:
        continue

    print(directory)
    species = species_names[directory.capitalize()]
    first = True

    # Loop through sam files in species directory
    for filename in os.listdir('mirgenedb_merged_extended_sam_filtered_no_outside_pos2_18_mismatch/' + str(directory)):
        if not filename.endswith('.sam'):
            continue
        print('\t' + filename)

    # Define filename as tissue (each original fasta file is one individual tissue)
        tissue = filename.split('/')[-1]

    # Parse through sam files, convert to pandas dataframe (only collapsed smallRNA seq samfiles, 
    # size is managable for pandas)
        print('\t parsing sam')
        df_ = parse_sam('mirgenedb_merged_extended_sam_filtered_no_outside_pos2_18_mismatch/' + str(directory) + '/' + str(filename), species, tissue)
        if first:
            df = df_
            first = False
        else:
            df = pd.concat([df, df_])

    # Write dataframe to .csv file 
    df = df.reset_index().drop('index', axis=1)
    print('writing ' + str(directory) + ' to .csv in: ' + str('mirgenedb_mismatch_df_directory/' + directory + '_mismatch_df.csv\n'))
    df.to_csv('mirgenedb_mismatch_df_directory/' + directory + '_mismatch_df.csv', sep='\t')



# Collapse the mismatch read files into one editing event / position per miRNA
print('\nCollapsing mismatch df file to only include one mismatch type_pos per miRNA')
for filename in os.listdir('mirgenedb_mismatch_df_directory/'):
    # Read file to pandas dataframe
    if not filename.endswith('.csv'):
        continue
    print(filename)
    df = pd.read_csv('mirgenedb_mismatch_df_directory/' + filename, sep='\t', index_col=0)

    # Sum up total number of reads per miRNA
    print('\tcount_total_reads')
    df['total'] = df['nr_reads'].groupby(df['miRNA']).transform('sum')

    # Transform dataframe to long format
    print('\ttransforming to long format')
    df = df.melt(id_vars=['ID', 'nr_reads', 'miRNA', 'arm', 'mature_star', 'species', 'tissue', 'read_start', 'read_seq', 'total'], 
                 value_vars=['from_to_1', 'from_to_2', 'from_to_3'], value_name='mismatch').drop('variable', axis=1)

    # Sum up total number for each mismatch per miRNA per position
    print('\tcounting_mismatches')
    total_mismatch = df.groupby(['miRNA', 'mismatch'])['nr_reads'].transform('sum').rename("total_mismatch").reset_index()
    df['total_mismatch'] = total_mismatch['total_mismatch']
    df = df[['miRNA', 'arm', 'mature_star', 'species', 'total', 'mismatch', 'total_mismatch']].drop_duplicates().dropna().sort_values(by=['total', 'total_mismatch'], ascending=False)

    # Edit the columns for easier downstream analysis
    print('\tsplit mismatch to-from and position into separate columns')
    df = df.reset_index().drop('index', axis=1)
    df['mismatch_pos'] = pd.Series(map(lambda x: split_mismatch_name(x)[1], df['mismatch']))
    df['edit_from']    = pd.Series(map(lambda x: split_edit_from_to(x)[0], df['mismatch']))
    df['edit_to']      = pd.Series(map(lambda x: split_edit_from_to(x)[1], df['mismatch']))
    df['mismatch']     = pd.Series(map(lambda x: split_mismatch_name(x)[0], df['mismatch']))
    df = df[['miRNA', 'arm', 'mature_star', 'species', 'total', 'mismatch', 'edit_from', 'edit_to', 'mismatch_pos', 'total_mismatch']]

    # Write the .csv file to specified directory
    print('\twriting ' + str('mirgenedb_collapsed_mismatch_df_directory/' + filename.split('_')[0] + '_collapsed_mismatch_df.csv'))
    df.to_csv('mirgenedb_collapsed_mismatch_df_directory/' + filename.split('_')[0] + '_collapsed_mismatch_df.csv', sep='\t')

# Collapse the mismatch read files into one editing event / position per miRNA per tissue
print('\nCollapsing mismatch df file to only include one mismatch type_pos per miRNA')
for filename in os.listdir('mirgenedb_mismatch_df_directory/'):
    # Read file to pandas dataframe
    if not filename.endswith('.csv'):
        continue
    print(filename)
    df = pd.read_csv('mirgenedb_mismatch_df_directory/' + filename, sep='\t', index_col=0)

    # Sum up total number of reads per miRNA
    print('\tcount_total_reads')
    df['total'] = df.groupby(['miRNA','tissue'])['nr_reads'].transform('sum')
    
    # Transform dataframe to long format
    print('\ttransforming to long format')
    df = df.melt(id_vars=['ID', 'nr_reads', 'miRNA', 'arm', 'mature_star', 'species', 'tissue', 'read_start', 'read_seq', 'total'],
                 value_vars=['from_to_1', 'from_to_2', 'from_to_3'], value_name='mismatch').drop('variable', axis=1)

    # Sum up total number for each mismatch per miRNA per position
    print('\tcounting_mismatches')
    total_mismatch = df.groupby(['miRNA', 'mismatch','tissue'])['nr_reads'].transform('sum').rename("total_mismatch").reset_index()
    df['total_mismatch'] = total_mismatch['total_mismatch']
    df = df[['miRNA', 'arm', 'mature_star', 'species', 'tissue', 'total', 'mismatch', 'total_mismatch']].drop_duplicates().dropna().sort_values(by=['total', 'total_mismatch'], ascending=False)

    # Edit the columns for easier downstream analysis
    print('\tsplit mismatch to-from and position into separate columns')
    df = df.reset_index().drop('index', axis=1)
    df['mismatch_pos'] = pd.Series(map(lambda x: split_mismatch_name(x)[1], df['mismatch']))
    df['edit_from']    = pd.Series(map(lambda x: split_edit_from_to(x)[0], df['mismatch']))
    df['edit_to']      = pd.Series(map(lambda x: split_edit_from_to(x)[1], df['mismatch']))
    df['mismatch']     = pd.Series(map(lambda x: split_mismatch_name(x)[0], df['mismatch']))
    df = df[['miRNA', 'arm', 'mature_star', 'species', 'tissue', 'total', 'mismatch', 'edit_from', 'edit_to', 'mismatch_pos', 'total_mismatch']]

    # Write the .csv file to specified directory
    print('\twriting ' + str('mirgenedb_collapsed_mismatch_by_tissue_df_directory/' + filename.split('_')[0] + '_collapsed_mismatch_by_tissue_df.csv'))
    df.to_csv('mirgenedb_collapsed_mismatch_by_tissue_df_directory/' + filename.split('_')[0] + '_collapsed_mismatch_by_tissue_df.csv', sep='\t')