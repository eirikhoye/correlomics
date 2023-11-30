import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
import seaborn as sns
import matplotlib.pyplot as plt
import os
from subprocess import *
from Bio import SeqIO
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import math
import copy
import json
import glob
import time

timestr = time.strftime("%Y%m%d-%H%M%S")

path = 'aligned_positions_extended_mature/'
li = []
print('loading data')
for directory in os.listdir(path):
    all_files = glob.glob(path + directory + "/*.csv")
    for filename in all_files:
        df = pd.read_csv(filename)
        li.append(df)
frame = pd.concat(li, axis=0, ignore_index=True)
frame.to_csv(path+'merged_aligned_positions_'+timestr+'.csv')

#bulges1 = pd.read_csv('results/bulges_copy.csv')
#bulges1 = bulges1.drop(['Unnamed: 0', 'Species'], axis=1)
#bulges2 = pd.read_csv('results/bulges_20210311-173602.csv')
#bulges2 = bulges2.drop(['Unnamed: 0', 'Species'], axis=1)
#bulges3 = pd.read_csv('results/bulges_20210312-092145.csv')
#bulges3 = bulges3.drop(['Unnamed: 0', 'Species'], axis=1)
#bulges4 = pd.concat([bulges1, bulges2, bulges3])
#bulges4.to_csv('results/test_bulges_4.csv')
bulges = pd.read_csv('results/bulges_20220621-135036.csv', sep=',')
bulges = bulges.drop(['Species', 'Unnamed: 0'], axis=1)
metadata = pd.read_csv('data/metadata_v7.csv', sep=';')
metadata = metadata.fillna(0)
frame_bulges = frame.merge(bulges, how='left', left_on='Gene_name', right_on='Gene_name')
timestr = time.strftime("%Y%m%d-%H%M%S")
frame_bulges_metadata = frame_bulges.merge(metadata, how='left', left_on='Gene_name', right_on='MirGeneDB_ID')
# DEBUG REMOVE
#frame_bulges_metadata.Mature_type = [n.lower() for n in frame_bulges_metadata.Mature_type]
# DEBUG REMOVE
frame_bulges_metadata.to_csv('results/frame_bulges_metadata_'+timestr+'.csv')
frame_bulges_metadata = frame_bulges_metadata.drop('Unnamed: 0', axis=1)

# A to I TEMP
path = 'mirgenedb_collapsed_mismatch_by_tissue_df_directory/'
is_first = True
for file in os.listdir(path):
    if not file.endswith('_collapsed_mismatch_by_tissue_df.csv'):
        continue
#    print(file)
    if is_first:
        mismatch_file = pd.read_csv(path+file, sep='\t')
        is_first = False
    else:
        mismatch_file_ = pd.read_csv(path+file, sep='\t')
        mismatch_file = pd.concat([mismatch_file, mismatch_file_], ignore_index=True).drop('Unnamed: 0', axis=1)

mismatch_subset = mismatch_file.loc[(mismatch_file.mismatch_pos.isin(list(range(2, 19))))]

mismatch_subset_atoi = mismatch_subset.loc[(mismatch_subset.mismatch == 'AtoG')]

### DEBUG REMOVE ###
mismatch_subset_atoi.to_csv('results/test_mismatch_subset_atoi.csv')
### DEBUG REMOVE ###

mismatch_subset_atoi['percentage'] = mismatch_subset_atoi.total_mismatch / mismatch_subset_atoi.total * 101
mismatch_subset_atoi['mean_percentage'] = mismatch_subset_atoi.groupby(by=['miRNA', 'tissue']).percentage.transform('mean')
mismatch_subset_atoi['max_percentage']  = mismatch_subset_atoi.groupby(by=['miRNA', 'tissue']).percentage.transform('max')
mismatch_subset_atoi['AtoI_max_pos'] = mismatch_subset_atoi.total_mismatch
mismatch_subset_atoi = mismatch_subset_atoi.drop('total_mismatch', axis=1)
mismatch_subset_atoi['total_AtoI'] = mismatch_subset_atoi.groupby(by=['miRNA', 'tissue']).AtoI_max_pos.transform('sum')
mismatch_subset_atoi = mismatch_subset_atoi[['miRNA', 'arm', 'mature_star', 'species', 'tissue', 'total', 'mismatch', 'edit_from', 'edit_to', 'mismatch_pos', 'percentage', 'mean_percentage', 'max_percentage', 'AtoI_max_pos', 'total_AtoI']]
mismatch_subset_atoi.columns = ['miRNA', 'arm', 'mature_star', 'species', 'tissue', 'total_filtered_reads', 'editing_event', 'edit_from', 'edit_to', 'mismatch_position', 'percentage_position', 'mean_percentage_miRNA', 'max_percentage_miRNA', 'AtoI_nr_reads_position', 'AtoI_nr_reads_total_miRNA']

mismatch_subset_atoi.to_csv('results/test_mismatch_subset_atoi.csv')

mismatch_subset_atoi_unique = mismatch_subset_atoi.loc[mismatch_subset_atoi.groupby(['miRNA', 'tissue'])['percentage_position'].idxmax()]
#mismatch_subset_atoi_unique = mismatch_subset_atoi.groupby(['miRNA', 'tissue']).apply(lambda group: group.nlargest(1, columns='percentage_position')).reset_index(level=-1, drop=True)
### DEBUG REMOVE ###
mismatch_subset_atoi_unique.to_csv('results/mismatch_subset_atoi_unique.csv')
### DEBUG REMOVE ###


# find the pesky comature...

mature = 'data/fasta/mature.fasta'
mature_dict = SeqIO.to_dict(SeqIO.parse(mature, 'fasta'))

comatures = []
patterns = [
'v2',
'v3',
'v4',
'v5',
'v6',
'v7',
'v8',
'v9'
]

for key in mature_dict.keys():
    mature_id = key.split('_')[0]
    if mature_id[-2:] in patterns:
        continue
    if mature_id + '_5p' and mature_id + '_3p' in mature_dict.keys():
        comatures.append(key)

for index, row in mismatch_subset_atoi_unique.iterrows():
    #print('print index, row', index, row)
    #print('print row.miRNA', row.miRNA)
    if row.miRNA in comatures and row.arm == '5p':
        mismatch_subset_atoi_unique.loc[index, 'mature_star'] = 'comature_5p'
    if row.miRNA in comatures and row.arm == '3p':
        mismatch_subset_atoi_unique.loc[index, 'mature_star'] = 'comature_3p'

### debug remove ###
mismatch_subset_atoi_unique.to_csv('results/mismatch_subset_atoi_unique_2.csv')
### debug_remove ###












mismatch_subset_atoi_unique_mature = mismatch_subset_atoi_unique[(mismatch_subset_atoi_unique.mature_star == 'mature') | (mismatch_subset_atoi_unique.mature_star == 'comature_5p') ]
mismatch_subset_atoi_unique_mature['miRNA'] = [i.split('_')[0] for i in mismatch_subset_atoi_unique_mature.miRNA]
mismatch_subset_atoi_unique_mature = mismatch_subset_atoi_unique_mature[['miRNA', 'species', 'tissue', 'total_filtered_reads', 'mismatch_position', 'AtoI_nr_reads_position', 'mean_percentage_miRNA', 'max_percentage_miRNA', 'AtoI_nr_reads_total_miRNA']]
mismatch_subset_atoi_unique_mature.columns = ['miRNA', 'species', 'tissue', 'total_mature_reads', 'mature_AtoI_max_pos', 'mature_max_position_nr_reads', 'mature_AtoI_mean_percentage', 'mature_AtoI_max_percentage', 'mature_AtoI_total_nr_reads']

### DEBUG REMOVE ###
mismatch_subset_atoi_unique_mature.to_csv('results/test_mismatch_subset_atoi_unique_mature.to_csv')
### DEBUG REMOVE ### 

mismatch_subset_atoi_unique_star = mismatch_subset_atoi_unique[(mismatch_subset_atoi_unique.mature_star == 'star') | (mismatch_subset_atoi_unique.mature_star == 'comature_3p')]
mismatch_subset_atoi_unique_star['miRNA'] = [i.split('_')[0] for i in mismatch_subset_atoi_unique_star.miRNA]
mismatch_subset_atoi_unique_star = mismatch_subset_atoi_unique_star[['miRNA', 'species', 'tissue', 'total_filtered_reads', 'mismatch_position', 'AtoI_nr_reads_position', 'mean_percentage_miRNA', 'max_percentage_miRNA', 'AtoI_nr_reads_total_miRNA']]
mismatch_subset_atoi_unique_star.columns = ['miRNA', 'species', 'tissue', 'total_star_reads', 'star_AtoI_max_pos', 'star_max_position_nr_reads', 'star_AtoI_mean_percentage', 'star_AtoI_max_percentage', 'star_AtoI_total_nr_reads']
mismatch_subset_atoi_unique_mature_star = mismatch_subset_atoi_unique_mature.merge(right=mismatch_subset_atoi_unique_star, how='outer', on=['miRNA', 'species', 'tissue']).fillna(0)
mismatch_subset_atoi_unique_mature_star.to_csv('../mirgenedb_editing_events/mirgenedb_collapsed_mismatch_by_tissue_df_directory/merged_species_collapsed_mismatch_df_AtoI_unique_mature_star.csv')
mismatch_subset_atoi_unique_mature_star.tissue = [i.split('.')[0] for i in mismatch_subset_atoi_unique_mature_star.tissue]
frame_bulges_metadata.Tissue = [i.split('.')[0] for i in frame_bulges_metadata.Tissue]
(frame_bulges_metadata)


### DEBUG REMOVE ###
mismatch_subset_atoi_unique_mature_star.to_csv('results/test_mismatch_subset_atoi_unique_mature_star.csv')
### DEBUG REMOVE ###

frame_bulges_metadata_atoi = frame_bulges_metadata.merge(mismatch_subset_atoi_unique_mature_star, how='left', left_on=['Gene_name', 'Tissue', 'Species'], right_on=['miRNA', 'tissue', 'species'])

### DEBUG REMOVE ###
frame_bulges_metadata_atoi.to_csv('results/test_frame_bulges_metadata_atoi.csv')
### DEBUG REMOVE ### 

frame_bulges_metadata_atoi = frame_bulges_metadata_atoi.drop(['miRNA', 'species', 'tissue'], axis=1).fillna(0)
dG_df = pd.read_csv('results/20220621-135036_deltaG.csv')
dG_df['Gene_name'] = [name.split('_')[0] for name in dG_df.miRNA]
dG_df = dG_df.drop('miRNA', axis=1)

frame_bulges_metadata_atoi_dG = frame_bulges_metadata_atoi.merge(dG_df, how='left', left_on='Gene_name', right_on='Gene_name')


stemfile    = open('data/fasta/stemfile.fasta', 'r')
stem        = SeqIO.to_dict(SeqIO.parse(stemfile, "fasta"))
prifile     = open('data/fasta/prifile.fasta', 'r')
pri         = SeqIO.to_dict(SeqIO.parse(prifile, "fasta"))

pre_df = pd.DataFrame({'Gene_name':[k.split('_')[0] for k in stem.keys()],
                       'pre_seq':  [str(v.seq).split('_')[0] for v in stem.values()]})

pri_df = pd.DataFrame({'Gene_name':[k.split('_')[0] for k in pri.keys()],
                       'pri_seq':  [str(v.seq).split('_')[0] for v in pri.values()]})

frame_bulges_metadata_atoi_dG = pd.merge(left=frame_bulges_metadata_atoi_dG, right=pre_df, left_on='Gene_name', right_on='Gene_name', how='left')
frame_bulges_metadata_atoi_dG = pd.merge(left=frame_bulges_metadata_atoi_dG, right=pri_df, left_on='Gene_name', right_on='Gene_name', how='left')


def count_len(sequence):
    if isinstance(sequence, str):
        return len(sequence)
    else:
        return("NA")



frame_bulges_metadata_atoi_dG['len_pre'] = list(map(count_len, frame_bulges_metadata_atoi_dG['pre_seq']))

frame_bulges_metadata_atoi_dG = frame_bulges_metadata_atoi_dG[['Gene_name',
 'Mature_Seq',
 'Star_Seq',
 'pre_seq',
 'len_pre',
 'pri_seq',
 'Mature_type',
 'Tissue',
 'Species',
 'Mature_Start_-5',
 'Mature_Start_-4',
 'Mature_Start_-3',
 'Mature_Start_-2',
 'Mature_Start_-1',
 'Mature_Start_0',
 'Mature_Start_1',
 'Mature_Start_2',
 'Mature_Start_3',
 'Mature_Start_4',
 'Mature_Start_5',
 'Mature_end_-5',
 'Mature_end_-4',
 'Mature_end_-3',
 'Mature_end_-2',
 'Mature_end_-1',
 'Mature_end_0',
 'Mature_end_1',
 'Mature_end_2',
 'Mature_end_3',
 'Mature_end_4',
 'Mature_end_5',
 'Star_Start_-5',
 'Star_Start_-4',
 'Star_Start_-3',
 'Star_Start_-2',
 'Star_Start_-1',
 'Star_Start_0',
 'Star_Start_1',
 'Star_Start_2',
 'Star_Start_3',
 'Star_Start_4',
 'Star_Start_5',
 'Star_End_-5',
 'Star_End_-4',
 'Star_End_-3',
 'Star_End_-2',
 'Star_End_-1',
 'Star_End_0',
 'Star_End_1',
 'Star_End_2',
 'Star_End_3',
 'Star_End_4',
 'Star_End_5',
 'Mature_total',
 'Star_total',
 '13.Length_loop',
 '14.Length_mature',
 '15.Length_star',
 '16.Number_of_bulges_mature',
 '17.Number_of_bulges_star',
 '18.Number_of_bulges_5p',
 '19.Number_of_bulges_3p',
 '20.Position_1_nucleotide_mature',
 '21.Position_1_nucleotide_star',
 '22.Position_1_nucleotide_5p',
 '23.Position_1_nucleotide_3p',
 '24.Position_1_bulge_mature',
 '25.Position_1_bulge_star',
 '26.Position_1_bulge_5p',
 '27.Position_1_bulge_3p',
 '28.SEED_mature',
 '29.SEED_star',
 '30.SEED_5p',
 '31.SEED_3p',
 'MirGeneDB_ID',
 'MiRBase_ID',
 'Family',
 'Seed',
 '5p_accession',
 '3p_accession',
 'Chromosome',
 'Start',
 'End',
 'Strand',
 'Node_of_origin_locus',
 'Node_of_origin_family',
 '3_NTU',
 ' UG ',
 'UGUG',
 'CNNC',
 'total_mature_reads',
 'mature_AtoI_max_pos',
 'mature_max_position_nr_reads',
 'mature_AtoI_mean_percentage',
 'mature_AtoI_max_percentage',
 'mature_AtoI_total_nr_reads',
 'total_star_reads',
 'star_AtoI_max_pos',
 'star_max_position_nr_reads',
 'star_AtoI_mean_percentage',
 'star_AtoI_max_percentage',
 'star_AtoI_total_nr_reads',
 'dG']]



print(timestr)
frame_bulges_metadata_atoi_dG.to_csv('results/frame_bulges_metadata_atoi_dG_'+timestr+'.csv')

for species in frame_bulges_metadata_atoi_dG.Species.unique():
    cm = frame_bulges_metadata_atoi_dG[frame_bulges_metadata_atoi_dG.Species == species][['Gene_name', 'Tissue', 'Mature_total']].pivot_table(
        index='Gene_name',
        columns='Tissue',
        values='Mature_total',
        aggfunc='first')
    cm.to_csv('results/count_matrixes/'+species+'_counts.csv')
    cm_rpm = cm/cm.sum()*1000000
    cm_rpm.to_csv('results/rpm_matrixes/'+species+'_rpms.csv')
































