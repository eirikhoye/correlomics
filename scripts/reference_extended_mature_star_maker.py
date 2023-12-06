from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import SeqIO
import os
import sys

"""
A python script to extend all MirGeneDB mature and star sequences
with 5nt 5prime and 3prime, from the precursor sequences.
"""

pri_mirna_dict = SeqIO.to_dict(SeqIO.parse(sys.argv[1], 'fasta'))
mature_star = sys.argv[2]
#logfile = open("/home/jcdenton/projects/mirgenedb_database/mirgenedb_merged_extended_genome/logfile.txt", 'w+') ### consider adding it back
outfile = open(sys.argv[3], 'w')

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

with open(mature_star, 'r') as file:
    for ID, SEQ in SimpleFastaParser(file):
        strand = ''
        if ID.endswith('*'):
            strand = 'star'
        if ID.endswith('p'):
            strand = 'mature'
        if ID.endswith('pre'):
            strand = 'pre'
        if strand == 'pre':
            continue
        if strand == 'mature' and ID[:-3] + '_pri' not in pri_mirna_dict:
            print(mature_star, ID, 'nope')
            #logfile.write(mature_star + ID + ' no_pri_sequence_found\n')
            continue
        if strand == 'star' and ID[:-4] + '_pri' not in pri_mirna_dict:
            print(mature_star, ID, 'nope')
            #logfile.write(mature_star + ID + ' no_pri_sequence_found\n')
            continue
        outfile.write('>' + ID + '\n')
        if strand == 'mature':
            start_pos = naive(str(SEQ), str(pri_mirna_dict[ID[:-3] + '_pri'].seq))['start'] - 5
            end_pos   = naive(str(SEQ), str(pri_mirna_dict[ID[:-3] + '_pri'].seq))['end'] + 5
            outfile.write(str(pri_mirna_dict[ID[:-3] + '_pri'].seq[start_pos : end_pos]) + '\n')
        if strand == 'star':
            start_pos = naive(str(SEQ), str(pri_mirna_dict[ID[:-4] + '_pri'].seq))['start'] - 5
            end_pos   = naive(str(SEQ), str(pri_mirna_dict[ID[:-4] + '_pri'].seq))['end'] + 5
            outfile.write(str(pri_mirna_dict[ID[:-4] + '_pri'].seq[start_pos : end_pos]) + '\n')


#for star in os.listdir('mirgenedb_star_genome/'):
#	outfile = open('mirgenedb_merged_extended_genome/' + str(star), 'a')
#	with open('mirgenedb_star_genome/' + star, 'r') as file:
#		for ID, SEQ in SimpleFastaParser(file):
#			if ID[:-4] + '_pri' not in pri_mirna_dict:
#				print('nope')
#				continue
#			outfile.write('>' + ID + '\n')
#			start_pos = naive(str(SEQ), str(pri_mirna_dict[ID[:-4] + '_pri'].seq))['start'] - 5
#			end_pos   = naive(str(SEQ), str(pri_mirna_dict[ID[:-4] + '_pri'].seq))['end'] + 5
#			outfile.write(str(pri_mirna_dict[ID[:-4] + '_pri'].seq[start_pos : end_pos]) + '\n')
