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
import time

def RNAfold(file_in):
    rnafold = Popen(["RNAfold", "-p", "-d2", "--noClosingGU", "--enforceConstraint", "-C"],
                    stdin=file_in,
                    stdout=PIPE,
                    universal_newlines=True)
    file_in.close()
    rnafold_output = rnafold.communicate()[0]
    output_lines = rnafold_output.splitlines()
    sequence = output_lines[0]
    structure = output_lines[1].split(None,1)[0].strip()
    energy = float(output_lines[1].rsplit("(",1)[1].strip("()"))
    return(sequence, structure, energy)


def naive(p, t):
    start = 0
    p = p.replace('U', 'T')
    t = t.replace('U', 'T')
    for i in range(len(t) - len(p) + 1):
        match = True
        for j in range(len(p)):
            if t[i+j] != p[j]:
                match = False
                break
        if match:
            start = i
    end = start + len(p)
    return({'start' : start, 'end' : end})

def CountWhobbles(input_fold):
    output_list = [0] * len(input_fold)
    for n in range(len(input_fold)):
        if n <= len(output_list):
            if input_fold[n] == '.':
                output_list[n] += 1
    return(output_list)

os.environ['PATH'] += os.pathsep + '/home/jcdenton/projects/mirgenedb_editing_events/src/RNAstructure/exe'

os.environ['DATAPATH'] = '/home/jcdenton/projects/mirgenedb_editing_events/src/RNAstructure/data_tables/'

##########################################################################################

# Load previous data

##########################################################################################

# Merge dataframe with MirGeneDB meta_data
#df = pd.read_csv('results/aligned_positions/merged_aligned_positions_20210423-175023.csv')

df = pd.read_csv('/home/jcdenton/projects/aligned_positions/merged_aligned_positions_extended_mature.csv')


#df.merge(metadata, left_on='Gene_name', right_on='MirGeneDB_ID', how='inner')
#df.to_csv('results/csv_file.csv')


##########################################################################################

# Create Bulge Data Frame

##########################################################################################

# load species name info
with open('species_dict.txt') as f:
    data = f.read()
species_names = json.loads(data)

# Define more column names
header_list = [
'Gene_name',
'Species',
'13.Length_loop',
'14.Length_mature',
'15.Length_star',
'16.Number_of_bulges_mature',
'17.Number_of_bulges_star',
'18.Number_of_bulges_5p',
'19.Number_of_bulges_3p',
'20.Position_1_nucleotide_mature', #(write A,C,G,T)
'21.Position_1_nucleotide_star',
'22.Position_1_nucleotide_5p',
'23.Position_1_nucleotide_3p',
'24.Position_1_bulge_mature', #(write 1 for bulge, 0 for no bulge)
'25.Position_1_bulge_star',
'26.Position_1_bulge_5p',
'27.Position_1_bulge_3p',
'28.SEED_mature', #(write_pos2-8:7nucleotides)
'29.SEED_star',
'30.SEED_5p',
'31.SEED_3p',
'32.Length_5p',
'33.Length_3p',
'34.mature_arm'
]



# Create dataframe for bulge data
bulges_df = pd.DataFrame(columns=header_list)
bulges_df['Gene_name'] = df['Gene_name'].unique()

stemfile = open('/home/jcdenton/projects/mirgenedb_database/merged_all-pri-30-30.fas', 'r')
loopfile = open('/home/jcdenton/projects/mirgenedb_database/loopfile.fas', 'r')
stem     = SeqIO.to_dict(SeqIO.parse(stemfile, "fasta"))
loop     = SeqIO.to_dict(SeqIO.parse(loopfile, "fasta"))

#prifile  = open('data/fasta/prifile.fasta', 'r')matures  = open('/home/jcdenton/projects/mirgenedb_database/mature.fas', 'r') 
#pri      = SeqIO.to_dict(SeqIO.parse(prifile, "fasta"))

deltaG_dict = {}

# define matures and co-matures
matures_5p = []
matures_3p = []
for n in matures:
    if n[-3:] == '5p\n':
        matures_5p.append(n[1:-4])
    if n[-3:] == '3p\n':
        matures_3p.append(n[1:-4])
comatures = []
for n in matures_5p:
    if n in matures_3p:
        comatures.append(n)
non_co_matures_5p = []
for n in matures_5p:
    if n not in comatures:
        non_co_matures_5p.append(n)
non_co_matures_3p = []
for n in matures_3p:
    if n not in comatures:
        non_co_matures_3p.append(n)


for n in stem:
    ### Loop over all pri-miRNA sequences in MirGeneDB
    gene_name = n[:-4]
    print(gene_name)
    #print(gene_name)
    species = species_names[stem[n].name[0:3]]

    ### DEBUG REMOVE 

#    if not species in ['Chicken']:
#        continue

    ### DEBUG REMOVE


    workingfile = open("workingfile.txt", "w")
    workingfile.write(str(stem[n].seq) + '\n')
    if stem[n].name[:-4] + '_loop' not in loop:
        continue
    loopstart = naive(str(loop[n[:-4] + '_loop'].seq), str(stem[n].seq))['start'] # find start position and end
    loopend   = naive(str(loop[n[:-4] + '_loop'].seq), str(stem[n].seq))['end']     # position of loop sequence
    workingfile.close()

    constraints = open("constraints.txt", "w")
    constraints.write('P' + ' ' + str(1)         + ' ' + str(0) + ' ' + str(17) + '\n')
#    constraints.write('P' + ' ' + str(loopstart + 1) + ' ' + str(0) + ' ' + str(loopend - loopstart - 2) + '\n')
    constraints.write('P' + ' ' + str(len(stem[n].seq) - 17 + 1) + ' ' + str(0) + ' ' + str(17))
    constraints.close()
    #workingfile = open("workingfile.txt", "r")
    os.system('mfold SEQ=workingfile.txt AUX=constraints.txt')
    os.system('ct2dot workingfile.txt.ct ALL workingfile.txt.dot')
    with open('workingfile.txt.dot', 'r') as file:
        dG = file.readline()
        seq = file.readline()
        dot = file.readline()
    deltaG_dict[stem[n].name] = dG
    folded = [seq.split('\n')[0], dot.split('\n')[0]]
    print(folded)
    file.close()
    i = CountWhobbles(folded[1])                    # Count '.' per position 
    UpBulgeLen   = len(i[30:loopstart])
    DownBulgeLen = len(i[loopend:-30])
    UpNrBulges   = sum(i[30:loopstart])
    DownNrBulges = sum(i[loopend:-30]) 

    #UpBulgeLenRatio = sum(i[30:loopstart]) / len(i[30:loopstart])
    #DownBulgeLenRatio = sum(i[loopend:-30]) / len(i[loopend:-30])
    UpSequence = stem[n].seq[:loopstart]
    DownSequence = stem[n].seq[loopend:]
    LoopSequence = stem[n].seq[loopstart:loopend]
    bulges_df.loc[ bulges_df.Gene_name == gene_name , '13.Length_loop'] = len(LoopSequence)
    bulges_df.loc[ bulges_df.Gene_name == gene_name, '18.Number_of_bulges_5p'] = UpNrBulges
    bulges_df.loc[ bulges_df.Gene_name == gene_name, '19.Number_of_bulges_3p'] = DownNrBulges
    bulges_df.loc[ bulges_df.Gene_name == gene_name, '32.Length_5p'] = UpBulgeLen
    bulges_df.loc[ bulges_df.Gene_name == gene_name, '33.Length_3p'] = DownBulgeLen
    if gene_name in matures_3p:
        bulges_df.loc[ bulges_df.Gene_name == gene_name, '14.Length_mature' ] = DownBulgeLen
        bulges_df.loc[ bulges_df.Gene_name == gene_name, '15.Length_star' ]   = UpBulgeLen
        bulges_df.loc[ bulges_df.Gene_name == gene_name, '16.Number_of_bulges_mature' ] = DownNrBulges
        bulges_df.loc[ bulges_df.Gene_name == gene_name, '17.Number_of_bulges_star' ]   = UpNrBulges
        bulges_df.loc[ bulges_df.Gene_name == gene_name, '20.Position_1_nucleotide_mature'] = str(stem[n].seq[loopend + 1]).replace('T', 'U')
        bulges_df.loc[ bulges_df.Gene_name == gene_name, '21.Position_1_nucleotide_star'] = str(stem[n].seq[30]).replace('T', 'U')
        bulges_df.loc[ bulges_df.Gene_name == gene_name, '34.mature_arm'] = 'matures_3p'
        if str(folded[1][loopend + 1]) == '.':
            bulges_df.loc[ bulges_df.Gene_name == gene_name, '24.Position_1_bulge_mature'] = 1
        else:
            bulges_df.loc[ bulges_df.Gene_name == gene_name, '24.Position_1_bulge_mature'] = 0
        if str(folded[1][30]) == '.':
            bulges_df.loc[ bulges_df.Gene_name == gene_name, '25.Position_1_bulge_star'] = 1
        else:
            bulges_df.loc[ bulges_df.Gene_name == gene_name, '25.Position_1_bulge_star'] = 0
        bulges_df.loc[ bulges_df.Gene_name == gene_name, '29.SEED_star' ] = str(stem[n].seq[31 : 38]).replace('T', 'U')
        bulges_df.loc[ bulges_df.Gene_name == gene_name, '28.SEED_mature' ] = str(stem[n].seq[loopend + 1 : loopend + 8]).replace('T', 'U')

    if gene_name in matures_5p:
        bulges_df.loc[ bulges_df.Gene_name == gene_name , '14.Length_mature'] = UpBulgeLen
        bulges_df.loc[ bulges_df.Gene_name == gene_name , '15.Length_star'] = DownBulgeLen
        bulges_df.loc[ bulges_df.Gene_name == gene_name, '16.Number_of_bulges_mature' ] = UpNrBulges
        bulges_df.loc[ bulges_df.Gene_name == gene_name, '17.Number_of_bulges_star' ]   = DownNrBulges
        bulges_df.loc[ bulges_df.Gene_name == gene_name, '20.Position_1_nucleotide_mature'] = str(stem[n].seq[30]).replace('T', 'U')
        bulges_df.loc[ bulges_df.Gene_name == gene_name, '21.Position_1_nucleotide_star'] = str(stem[n].seq[loopend]).replace('T', 'U')
        bulges_df.loc[ bulges_df.Gene_name == gene_name, '34.mature_arm'] = 'matures_5p'

        if str(folded[1][30]) == '.':
            bulges_df.loc[ bulges_df.Gene_name == gene_name, '24.Position_1_bulge_mature'] = 1
        else:
            bulges_df.loc[ bulges_df.Gene_name == gene_name, '24.Position_1_bulge_mature'] = 0
        if str(folded[1][loopend + 1]) == '.':
            bulges_df.loc[ bulges_df.Gene_name == gene_name, '25.Position_1_bulge_star'] = 1
        else:
            bulges_df.loc[ bulges_df.Gene_name == gene_name, '25.Position_1_bulge_star'] = 0
        bulges_df.loc[ bulges_df.Gene_name == gene_name, '28.SEED_mature' ] = str(stem[n].seq[31 : 38]).replace('T', 'U')
        bulges_df.loc[ bulges_df.Gene_name == gene_name, '29.SEED_star' ] = str(stem[n].seq[loopend + 1 : loopend + 8]).replace('T', 'U')

    bulges_df.loc[ bulges_df.Gene_name == gene_name, '18.Number_of_bulges_5p' ] = UpNrBulges
    bulges_df.loc[ bulges_df.Gene_name == gene_name, '19.Number_of_bulges_3p' ] = DownNrBulges
    bulges_df.loc[ bulges_df.Gene_name == gene_name, '32.Length_5p'] = UpBulgeLen
    bulges_df.loc[ bulges_df.Gene_name == gene_name, '33.Length_3p'] = DownBulgeLen

    bulges_df.loc[ bulges_df.Gene_name == gene_name, '22.Position_1_nucleotide_5p'] = str(stem[n].seq[30]).replace('T', 'U')
    bulges_df.loc[ bulges_df.Gene_name == gene_name, '23.Position_1_nucleotide_3p'] = str(stem[n].seq[loopend]).replace('T', 'U')
    if str(folded[1][30]) == '.':
        bulges_df.loc[ bulges_df.Gene_name == gene_name, '26.Position_1_bulge_5p'] = 1
    else:
        bulges_df.loc[ bulges_df.Gene_name == gene_name, '26.Position_1_bulge_5p'] = 0
    if str(folded[1][loopend + 1]) == '.':
        bulges_df.loc[ bulges_df.Gene_name == gene_name, '27.Position_1_bulge_3p'] = 1
    else:
        bulges_df.loc[ bulges_df.Gene_name == gene_name, '27.Position_1_bulge_3p'] = 0
    bulges_df.loc[ bulges_df.Gene_name == gene_name, '30.SEED_5p' ] = str(stem[n].seq[31 : 38]).replace('T', 'U')
    bulges_df.loc[ bulges_df.Gene_name == gene_name, '31.SEED_3p' ] = str(stem[n].seq[loopend + 1 : loopend + 8]).replace('T', 'U')
    
timestr = time.strftime("%Y%m%d-%H%M%S")
#print(bulges_df.head())
bulges_df.to_csv('results/bulges_'+ timestr +'.csv')
dG_df = pd.DataFrame.from_dict({key:value[6:12] for key,value in deltaG_dict.items()}, orient='index').reset_index()
dG_df.columns = ['miRNA', 'dG']
dG_df.to_csv('results/'+timestr+'_deltaG.csv', index=False)



