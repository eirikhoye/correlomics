import pandas as pd
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
import matplotlib.pyplot as plt
import os
from subprocess import *
from Bio import SeqIO
import pickle

stemfile   = open('/home/jcdenton/projects/mirgenedb_database/merged_all-pri-30-30.fas', 'r')
maturefile = open('/home/jcdenton/projects/mirgenedb_database/mature.fas', 'r') 
starfile   = open('/home/jcdenton/projects/mirgenedb_database/star.fas')
stem   = SeqIO.to_dict(SeqIO.parse(stemfile, "fasta"))
mature = SeqIO.to_dict(SeqIO.parse(maturefile, "fasta"))
star   = SeqIO.to_dict(SeqIO.parse(starfile, "fasta"))

os.environ['PATH'] += os.pathsep + '/home/jcdenton/projects/mirgenedb_editing_events/src/RNAstructure/exe'
os.environ['DATAPATH'] = '/home/jcdenton/projects/mirgenedb_editing_events/src/RNAstructure/data_tables/'

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

class foldings():
    def __init__(self, mirgenedb_id):
        self.mirgenedb_id = mirgenedb_id
        self.stem, self.mature, self.star, self._5p_arm, self._3p_arm = self.fetch_sequences()
        
    
    
    def is_5p_or_3p_mature(self):
        mature_id = [id for id in mature]
        if self.mirgenedb_id + "_5p" in mature_id and self.mirgenedb_id + "_3p" in mature_id:
            return "co_mature"
        elif self.mirgenedb_id + "_5p" in mature_id:
            return "5p_mature"    
        else:
            return "3p_mature"
        

    def fetch_sequences(self):
        mature_arm = self.is_5p_or_3p_mature()
        if mature_arm == "5p_mature":
            return (str(stem[self.mirgenedb_id + "_pri"].seq), 
                    str(mature[self.mirgenedb_id + "_5p"].seq), 
                    str(star[self.mirgenedb_id + "_3p*"].seq),
                    str(mature[self.mirgenedb_id +"_5p"].seq),
                    str(star[self.mirgenedb_id + "_3p*"].seq))
        if mature_arm == "3p_mature":
            return (str(stem[self.mirgenedb_id + "_pri"].seq), 
                    str(mature[self.mirgenedb_id + "_3p"].seq), 
                    str(star[self.mirgenedb_id + "_5p*"].seq),
                    str(star[self.mirgenedb_id + "_5p*"].seq),
                    str(mature[self.mirgenedb_id + "_3p"].seq))
        if mature_arm == 'co_mature':
            return (str(stem[self.mirgenedb_id + "_pri"].seq), 
                    str(mature[self.mirgenedb_id + "_5p"].seq), 
                    str(mature[self.mirgenedb_id + "_3p"].seq),
                    str(mature[self.mirgenedb_id + "_5p"].seq),
                    str(mature[self.mirgenedb_id + "_3p"].seq))


    def fold_hairpin(self):

        self.end_5p_arm = naive(self._5p_arm, self.stem)['end']   # get mature end position
        self.start_3p_arm = naive(self._3p_arm, self.stem)['start']   # get position of first star
        
        os.chdir('tmp')
        workingfile = open('workingfile.txt', 'w')
        workingfile.write(self.stem + '\n')
        workingfile.close()
        constraints = open("constraints.txt", "w")
        constraints.write('P' + ' ' + str(1)         + ' ' + str(0) + ' ' + str(17) + '\n')     # define 5p ds
#       constraints.write('P' + ' ' + str(loopstart + 1) + ' ' + str(0) + ' ' + str(loopend - loopstart - 2) + '\n') # define loop ds
        constraints.write('P' + ' ' + str(len(self.stem) - 17 + 1) + ' ' + str(0) + ' ' + str(17))    # define 3p ds
        constraints.close()

        os.system('mfold SEQ=workingfile.txt AUX=constraints.txt')
        os.system('ct2dot workingfile.txt.ct ALL workingfile.txt.dot')
        with open('workingfile.txt.dot', 'r') as file:
            self.dG = file.readline()
            self.seq = file.readline()
            self.dot = file.readline()

        os.chdir('..')

    def find_motif(self):
        self.start_5p_arm = naive(self._5p_arm, self.stem)['start']   # get mature end position
        self.end_3p_arm = naive(self._3p_arm, self.stem)['end']   # get position of first star

        pos_5p = [i + self.start_5p_arm for i in [-6,-5, -4, -3, -2, -1]]
        pos_3p = [i + self.end_3p_arm for i in [+3, +2, +1, 0]]

        dicer_cleave = [
            [self.hairpin['seq'][i] for i in pos_5p],
            [self.hairpin['dot'][i] for i in pos_5p],
            [self.hairpin['dot'][i] for i in pos_3p],
            [self.hairpin['seq'][i] for i in pos_3p]
                        ]


        self.GHG = [0,0,0]
        if dicer_cleave[0][1] == 'C':
            self.GHG[0] = 1
        if dicer_cleave[3][1] == 'G':
            self.GHG[1] = 1
        if dicer_cleave[1][0] == '.' and dicer_cleave[2][0] == '.':
            self.GHG[2] = 1

        

    def fetch_hairpin(self):
        with open(hairpin_dict + self.mirgenedb_id + '.pkl', 'rb') as f:
            self.hairpin = pickle.load(f)


hairpin_dict = '/home/jcdenton/projects/mirgenedb_database/hairpins_restrict_loop/'
#### DEBUG ####
hairpin_dict = '/home/jcdenton/projects/correlomics/test_folds/'
#### DEBUG ####

"""
test_miR = "Hsa-Mir-26-P1"

test = foldings(test_miR)
test.fetch_hairpin()
test.find_motif()
GHG_dict = {}

GHG_dict[test_miR] = test.GHG
print(GHG_dict)
"""



GHG_dict = {}

for file in os.listdir(hairpin_dict):
    mirgenedb_id = file[:-4]
    fold = foldings(mirgenedb_id)
    fold.fetch_hairpin()
    fold.find_motif()
    GHG_dict[mirgenedb_id] = fold.GHG


df = pd.DataFrame.from_dict(GHG_dict, orient='index')
df.columns = ['C', 'G', 'm']
df.to_csv('mGHG_motifs.csv')
