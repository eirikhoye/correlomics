import sys
import re
import os
import numbers
import numpy as np

"""
Python script, takes sam files as input, outputs sam files were read
alignments are filtered. 

Filtering is done to prioritise reads aligning to multiple miRNA
reference mature/star sequences.

Filtering is score based, with the following criteria:
    Reads that did not align:                 Score = 0
    Reads that did align to reference:        Score = 100
    If alignment has mismatch in pos 2-19 nt: Score - 999
    If alignment has 1 mismatch:              Score - 10
    If alignment has 2 mismatch:              Score - 20
    If alignment has 3 mismatch:              Score - 30

Filtering keeps only the read or reads with the highest score.

Outputs to 
mirgenedb_merged_extended_sam_filtered_no_outside_pos2_18_mismatch
"""

input_file = sys.argv[1]
output_file = sys.argv[2]

readID = {}
refReads = {}
with open(input_file, 'r') as file:

	for line in file:
		if '@HD' in line:
			header = line
			continue
		if 'SQ' in line:
			refReads[line] = line
			continue
		if '@PG' in line:
			parameters = line
			continue
		line = line.split('\t')
		if line[2] == '*':
			continue
		if line[0] not in readID:
			readID[line[0]] = []
			readID[line[0]].append(line)
		else:
			readID[line[0]].append(line)

for key in readID:
	score ={}
	for element in readID[key]:
		if element[2] == '*':
			element[11] = element[11].rstrip()
			score[str(element)] = 0
			continue
		
		nr_mismatch = int(element[13].split(':')[2].split("\n")[0])
		sequence = element[9]
		read_start = int(element[3])
		mismatch_pos = element[12].split(':')[2]
		mismatch_pos = re.split('(\d+)', mismatch_pos)
		k = [1, 3, 5]
		y = [0] * nr_mismatch
		for i in range(nr_mismatch):
			if i > 0:
				y[i] = read_start - 5 + int(mismatch_pos[k[i]]) + int(y[i-1])
			else:
				y[0] = read_start - 5 + int(mismatch_pos[k[i]])

		element[13] = element[13].rstrip()


		score[str(element)] = 100
		seed = list(range(2, 19))

		for i in y:
			if i not in seed:
				score[str(element)] -= 999
		if element[13] == 'NM:i:0':
			score[str(element)] += 0

		if element[13] == 'NM:i:1':
			score[str(element)] -= 10

		if element[13] == 'NM:i:2':
			score[str(element)] -= 20

		if element[13] == 'NM:i:3':
			score[str(element)] -= 30

	mismatch_in_seed = False
	def keywithmaxval(d):
		# Return list of reads with best score
		highest = max(d.values())
		if highest < 0:
			global mismatch_in_seed
			mismatch_in_seed = True
		return([k for k, v in d.items() if v == highest])
	topread = keywithmaxval(score)
	topread = "\t".join(topread)
	for ch in ['[',']',"'", " "]:
			if ch in topread:
				topread = topread.replace(ch,'')
	topread = topread.rstrip().replace('\t', '\n').replace(',','\t')
	readID[key] = topread
	if mismatch_in_seed == True:
		readID[key] = ''

with open(output_file, 'w') as outfile:
	outfile.write(header)
	for key in refReads:
		outfile.write(refReads[key])
	outfile.write(parameters)
	for key in readID:
		if readID[key] == '':
			continue
		outfile.write(readID[key] + '\n')
