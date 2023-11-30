from subprocess import *
from Bio import SeqIO
import os

os.environ['PATH'] += os.pathsep + '/home/jcdenton/projects/mirgenedb_editing_events/src/RNAstructure/exe/'
os.environ['DATAPATH'] = '/home/jcdenton/projects/mirgenedb_editing_events/src/RNAstructure/data_tables/'
print(os.environ['PATH'])


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



class foldings():
    def __init__(self, mirgenedb_id):
        self.mirgenedb_id = mirgenedb_id
        #self.pri_mirna = # make function to fetch primirna

    #def get_sequences():
        
    
    #def fold():

    #    workingfile = open("tmp/workingfile.txt", "w")
    #    workingfile.write(str(self.pri_mirna))
