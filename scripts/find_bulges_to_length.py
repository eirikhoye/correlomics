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
        self.end_5p_arm = naive(self._5p_arm, self.stem)['end']   # get mature end position
        self.start_3p_arm = naive(self._3p_arm, self.stem)['start']   # get position of first star

        pos_5p = [i + self.end_5p_arm for i in [-4, -3, -2, -1]]
        pos_3p = [i + self.start_3p_arm for i in [1, 0, -1, -2]]

        dicer_cleave = [
            [self.hairpin['seq'][i] for i in pos_5p],
            [self.hairpin['dot'][i] for i in pos_5p],
            [self.hairpin['dot'][i] for i in pos_3p],
            [self.hairpin['seq'][i] for i in pos_3p]
                        ]

        self.GYM = [0,0,0]
        if dicer_cleave[3][-2] == 'G':
            self.GYM[0] = 1
        if dicer_cleave[3][-3] == 'C' or dicer_cleave[3][-3] == 'T' or dicer_cleave[3][-3] == 'U':
            self.GYM[1] = 1
        if dicer_cleave[2][-4] == '.':
            self.GYM[2] = 1
    
    def count_bulges(self, fold_string):
        nr_bulges = 0
        for i in fold_string:
            if i == '.':
                nr_bulges += 1
        return(nr_bulges)

    def count_basepair(self, fold_string):
        nr_basepair = 0
        for i in fold_string:
            if i == '(' or i == ')':
                nr_basepair += 1
        return(nr_basepair)

    def check_misfold(self):
        self.end_5p_arm = naive(self._5p_arm, self.stem)['end']   # get mature end position
        self.start_3p_arm = naive(self._3p_arm, self.stem)['start']   # get position of first star
        misfold = False
        for i in range(len(self.hairpin['dot'])):
            if (str(self.hairpin['dot'][i]) == ')' and i < self.end_5p_arm):
                misfold = True
            elif (str(self.hairpin['dot'][i]) == '(' and i > self.start_3p_arm):
                misfold = True

        return misfold

    def find_bulge_length(self):
        pos_5p_arm = naive(self._5p_arm, self.stem)
        pos_3p_arm = naive(self._3p_arm, self.stem)
        
        bulges_5p = self.hairpin['dot'][pos_5p_arm['start']:pos_5p_arm['end']]
        bulges_3p = self.hairpin['dot'][pos_3p_arm['start']:pos_3p_arm['end']-2] # Get correct counting

        self.bulge_len_dict = {
                '5p_length' : len(self._5p_arm),
                '5p_bulges' : self.count_bulges(bulges_5p),
                '5p_bulges_len' : float(self.count_bulges(bulges_5p)) / (len(self._5p_arm)),
                '5p_basepair' : self.count_basepair(bulges_5p),
                '5p_basepair_len' : float(self.count_basepair(bulges_5p)) / (len(self._5p_arm)),
                '3p_length' : len(self._3p_arm),
                '3p_bulges' : self.count_bulges(bulges_3p),
                '3p_bulges_len' : float(self.count_bulges(bulges_3p)) / len(self._3p_arm),
                '3p_basepair' : self.count_basepair(bulges_3p),
                '3p_basepair_len' : float(self.count_basepair(bulges_3p)) / (len(self._3p_arm)) 
                }

    def fetch_hairpin(self):
        with open(hairpin_dict + self.mirgenedb_id + '.pkl', 'rb') as f:
            self.hairpin = pickle.load(f)


hairpin_dict = '/home/jcdenton/projects/mirgenedb_database/hairpins_PRE_no_constraints/'
#### DEBUG ####
hairpin_dict = '/home/jcdenton/projects/correlomics/test_folds/'
#### DEBUG ####


"""
test_miR = "Hsa-Mir-3174"

test = foldings(test_miR)
test.fetch_hairpin()
#test.find_bulge_length()
#GYM_dict = {}

#GYM_dict[test_miR] = test.GYM
#print(GYM_dict)
#print(test.hairpin)
#print(test.hairpin['dot'])
#print(test.hairpin['dot'][30])
##print(test.hairpin['dot'][:-30])
print(test.find_bulge_length())
"""

bulge_len_dict = {}
misfold_file = open('misfolds.txt', 'w')

for file in os.listdir(hairpin_dict):
    mirgenedb_id = file[:-4]
    fold = foldings(mirgenedb_id)
    fold.fetch_hairpin()
    if fold.check_misfold() == True:
        misfold_file.write(mirgenedb_id+'\n')
        misfold_file.write(fold.hairpin['seq'])
        misfold_file.write(fold.hairpin['dot'])
        continue
    fold.find_bulge_length()
    bulge_len_dict[mirgenedb_id] = fold.bulge_len_dict
misfold_file.close()

df = pd.DataFrame.from_dict(bulge_len_dict, orient='index')
df.columns = ['5p_length', '5p_bulges', '5p_bulges_len', '5p_basepair', '5p_basepair_len',
              '3p_length', '3p_bulges', '3p_bulges_len', '3p_basepair', '3p_basepair_len'
             ]
df.to_csv('bulge_len_ratio.csv')




"""



hairpin_dict = '/home/jcdenton/projects/mirgenedb_database/hairpins_PRE/'

test_miR = "Hsa-Mir-3174"

test = foldings(test_miR)
test.fetch_hairpin()
#test.find_bulge_length()
#GYM_dict = {}

#GYM_dict[test_miR] = test.GYM
#print(GYM_dict)
#print(test.hairpin)
#print(test.hairpin['dot'])
#print(test.hairpin['dot'][30])
##print(test.hairpin['dot'][:-30])
print(test.find_bulge_length())


bulge_len_dict = {}
misfold_file = open('misfolds.txt', 'w')

for file in os.listdir(hairpin_dict):
    mirgenedb_id = file[:-4]
    fold = foldings(mirgenedb_id)
    fold.fetch_hairpin()
    if fold.check_misfold() == True:
        misfold_file.write(mirgenedb_id+'\n')
        misfold_file.write(fold.hairpin['seq'])
        misfold_file.write(fold.hairpin['dot'])
        continue
    fold.find_bulge_length()
    bulge_len_dict[mirgenedb_id] = fold.bulge_len_dict
misfold_file.close()

df = pd.DataFrame.from_dict(bulge_len_dict, orient='index')
df.columns = ['5p_length', '5p_bulges', '5p_bulges_len', '5p_basepair', '5p_basepair_len'
#              '3p_length', '3p_bulges', '3p_bulges_len'
             ]
df.to_csv('bulge_len_ratio_restricted_loop_PRE.csv')

"""