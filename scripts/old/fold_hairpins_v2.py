import pandas as pd
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
import matplotlib.pyplot as plt
import os
from subprocess import *
from Bio import SeqIO
import pickle
import argparse

"""
Define input parameters
"""
parser = argparse.ArgumentParser()
parser.add_argument("--verbosity", help="increase output verbosity")
args = parser.parse_args()
if args.verbosity:
    print('verbosity turned on')

"""
Define input files
"""
stemfile   = open('/home/jcdenton/projects/mirgenedb_database/merged_all-pri-30-30.fas', 'r')
maturefile = open('/home/jcdenton/projects/mirgenedb_database/mature.fas', 'r') 
starfile   = open('/home/jcdenton/projects/mirgenedb_database/star.fas')
loopfile   = open('/home/jcdenton/projects/mirgenedb_database/loopfile.fas')
stem   = SeqIO.to_dict(SeqIO.parse(stemfile, "fasta"))
mature = SeqIO.to_dict(SeqIO.parse(maturefile, "fasta"))
star   = SeqIO.to_dict(SeqIO.parse(starfile, "fasta"))
loop   = SeqIO.to_dict(SeqIO.parse(loopfile, "fasta"))
outdir = '/home/jcdenton/projects/mirgenedb_database/hairpins_pri/'

"""
define path variables for mfold
"""
os.environ['PATH'] += os.pathsep + '/home/jcdenton/projects/mirgenedb_editing_events/src/RNAstructure/exe'
os.environ['DATAPATH'] = '/home/jcdenton/projects/mirgenedb_editing_events/src/RNAstructure/data_tables/'

def naive(p, t):
    """
    Get start and end position for query string p in target string t
    """
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
        
    def fetch_loop(self):
        if self.mirgenedb_id + '_loop' in loop:
            self.loop = str(loop[self.mirgenedb_id + '_loop'].seq)
        
    def fold_hairpin(self, constrain_stem = False, constrain_loop = False):
        """
        TODO
        add function to check whether mirna has loop sequence, return error if not
        
        """
        self.fetch_loop()
        self.end_5p_arm   = naive(self._5p_arm, self.stem)['end']   # get mature end position
        self.start_3p_arm = naive(self._3p_arm, self.stem)['start']   # get position of first star
        self.start_loop   = naive(self.loop, self.stem)['start']
        self.end_loop     = naive(self.loop, self.stem)['end']
        os.chdir('tmp')
        workingfile = open('workingfile.txt', 'w')
        workingfile.write(self.stem + '\n')
        workingfile.close()
        constraints = open("constraints.txt", "w")
        if constrain_stem:
            constraints.write('P' + ' ' + str(1) + ' ' + str(0) + ' ' + str(17) + '\n')     # define 5p ds
        if constrain_loop:
            constraints.write('P' + ' ' + str(self.start_loop + 1) + ' ' + str(0) + ' ' + str(self.end_loop - self.start_loop - 2) + '\n') 
        if constrain_stem:
            constraints.write('P' + ' ' + str(len(self.stem) - 17 + 1) + ' ' + str(0) + ' ' + str(17))    # define 3p ds
        constraints.close()

        os.system('mfold SEQ=workingfile.txt AUX=constraints.txt')
        os.system('ct2dot workingfile.txt.ct ALL workingfile.txt.dot')
        with open('workingfile.txt.dot', 'r') as file:
            self.dG = file.readline()
            self.seq = file.readline()
            self.dot = file.readline()

        os.chdir('..')

        self.fold_dict = {'dG' : self.dG, 'seq' : self.seq, 'dot' : self.dot}

    def check_misfold(self, stem_fold): # Check if there are opposite folds before or after center of loop
        misfold = False
        loop_middle = int(self.start_3p_arm - ((self.start_3p_arm - self.end_5p_arm) / 2))
        for i in range(len(stem_fold)):
            if (str(stem_fold[i]) == ')' and i < loop_middle):
                misfold = True
            elif (str(stem_fold[i]) == '(' and i > loop_middle):
                misfold = True
        return misfold

misfold_no_constraints_file = open('misfold_no_constraints.txt', 'w')
misfold_constrain_lower_stem_file = open('misfold_constrain_lower_stem.txt', 'w')
misfold_constrain_lower_stem_and_loop_file = open('misfold_constrain_lower_stem_and_loop.txt', 'w')

for id in stem.keys():
    if "Mir-451" in id:
        continue
    hairpin = foldings(id[:-4]) # drop miRNA without loop seq
    if hairpin.mirgenedb_id + '_loop' not in loop.keys():
        continue
    hairpin.fold_hairpin()  # no constraints
    if hairpin.check_misfold(hairpin.fold_dict['dot']) == False:
        with open(outdir + hairpin.mirgenedb_id + '.pkl', 'wb') as f:
            pickle.dump(hairpin.fold_dict, f)
    else:
        misfold_no_constraints_file.write(id+'\n')
        misfold_no_constraints_file.write(hairpin.fold_dict['seq'])
        misfold_no_constraints_file.write(hairpin.fold_dict['dot'])
        hairpin.fold_hairpin(constrain_stem = True) # Constrain lower stem to not fold
        if hairpin.check_misfold(hairpin.fold_dict['dot']) == False:
            with open(outdir + hairpin.mirgenedb_id + '.pkl', 'wb') as f:
                pickle.dump(hairpin.fold_dict, f)
        else:
            misfold_constrain_lower_stem_file.write(id+'\n')
            misfold_constrain_lower_stem_file.write(hairpin.fold_dict['seq'])
            misfold_constrain_lower_stem_file.write(hairpin.fold_dict['dot'])
            hairpin.fold_hairpin(constrain_stem = True, constrain_loop = True)  # Constrain lower stem and loop to not fold

            if hairpin.check_misfold(hairpin.fold_dict['dot']) == False:
                with open(outdir + hairpin.mirgenedb_id + '.pkl', 'wb') as f:
                    pickle.dump(hairpin.fold_dict, f)
            else:
                misfold_constrain_lower_stem_and_loop_file.write(id+'\n')
                misfold_constrain_lower_stem_and_loop_file.write(hairpin.fold_dict['seq'])
                misfold_constrain_lower_stem_and_loop_file.write(hairpin.fold_dict['dot'])
                
misfold_no_constraints_file.close()
misfold_constrain_lower_stem_file.close()
misfold_constrain_lower_stem_and_loop_file.close()