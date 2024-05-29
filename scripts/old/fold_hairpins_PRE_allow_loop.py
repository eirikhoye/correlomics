import pandas as pd
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
import matplotlib.pyplot as plt
import os
from subprocess import *
from Bio import SeqIO
import pickle

stemfile   = open('/home/jcdenton/projects/mirgenedb_database/precursor.fas', 'r')
maturefile = open('/home/jcdenton/projects/mirgenedb_database/mature.fas', 'r') 
starfile   = open('/home/jcdenton/projects/mirgenedb_database/star.fas')
loopfile   = open('/home/jcdenton/projects/mirgenedb_database/loopfile.fas')
stem   = SeqIO.to_dict(SeqIO.parse(stemfile, "fasta"))
mature = SeqIO.to_dict(SeqIO.parse(maturefile, "fasta"))
star   = SeqIO.to_dict(SeqIO.parse(starfile, "fasta"))
loop   = SeqIO.to_dict(SeqIO.parse(loopfile, "fasta"))

os.environ['PATH'] += os.pathsep + '/home/jcdenton/projects/mirgenedb_editing_events/src/RNAstructure/exe'
os.environ['DATAPATH'] = '/home/jcdenton/projects/mirgenedb_editing_events/src/RNAstructure/data_tables/'

outdir = '/home/jcdenton/projects/mirgenedb_database/hairpins_PRE_no_constraints/'

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
            return (str(stem[self.mirgenedb_id + "_pre"].seq), 
                    str(mature[self.mirgenedb_id + "_5p"].seq), 
                    str(star[self.mirgenedb_id + "_3p*"].seq),
                    str(mature[self.mirgenedb_id +"_5p"].seq),
                    str(star[self.mirgenedb_id + "_3p*"].seq))
        if mature_arm == "3p_mature":
            return (str(stem[self.mirgenedb_id + "_pre"].seq), 
                    str(mature[self.mirgenedb_id + "_3p"].seq), 
                    str(star[self.mirgenedb_id + "_5p*"].seq),
                    str(mature[self.mirgenedb_id + "_3p"].seq),
                    str(star[self.mirgenedb_id + "_5p*"].seq))
        if mature_arm == 'co_mature':
            return (str(stem[self.mirgenedb_id + "_pre"].seq), 
                    str(mature[self.mirgenedb_id + "_5p"].seq), 
                    str(mature[self.mirgenedb_id + "_3p"].seq),
                    str(mature[self.mirgenedb_id + "_5p"].seq),
                    str(mature[self.mirgenedb_id + "_3p"].seq))
        
    def fetch_loop(self):
        if self.mirgenedb_id + '_loop' in loop:
            self.loop = str(loop[self.mirgenedb_id + '_loop'].seq)
        


    def fold_hairpin(self):
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
        #constraints.write('P' + ' ' + str(1)         + ' ' + str(0) + ' ' + str(17) + '\n')     # define 5p ds
#       constraints.write('P' + ' ' + str(loopstart + 1) + ' ' + str(0) + ' ' + str(loopend - loopstart - 2) + '\n') # define loop ds
        #constraints.write('P' + ' ' + str(self.start_loop + 1) + ' ' + str(0) + ' ' + str(self.end_loop - self.start_loop - 2) + '\n') 
        #constraints.write('P' + ' ' + str(len(self.stem) - 17 + 1) + ' ' + str(0) + ' ' + str(17))    # define 3p ds
        constraints.close()

        os.system('mfold SEQ=workingfile.txt AUX=constraints.txt')
        os.system('ct2dot workingfile.txt.ct ALL workingfile.txt.dot')
        with open('workingfile.txt.dot', 'r') as file:
            self.dG = file.readline()
            self.seq = file.readline()
            self.dot = file.readline()

        os.chdir('..')

        self.fold_dict = {'dG' : self.dG, 'seq' : self.seq, 'dot' : self.dot}
        with open(outdir + self.mirgenedb_id + '.pkl', 'wb') as f:
            pickle.dump(self.fold_dict, f)

for id in stem.keys():
    if "Mir-451" in id:
        continue
    hairpin = foldings(id[:-4])                     # drop miRNA without loop seq
    if hairpin.mirgenedb_id + '_loop' not in loop.keys():
        continue
    hairpin.fold_hairpin()
