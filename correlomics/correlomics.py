import pandas as pd
from subprocess import *
from Bio import SeqIO
import json
import glob
import os

class Correlomics():
	def __init__(self, species):
		"""
		Initialize the class with the species name,
		and parse the reference files
		"""
		self.species = species
		species_names = self.load_species_dict()
		self.reference('data/fasta/prifile.fasta',
						'data/fasta/mature_extended.fas',
						'data/fasta/star_extended_fake_451.fas',
						'data/fasta/mature.fasta',
						'data/fasta/starfile_fake_451.fasta', self.species)


	def load_species_dict(self):
		"""Load species dictionary from file"""
		with open('species_dict.txt') as f:
			data = f.read()
		return json.loads(data)


	def reference(self, pri_ref, mature_ext, star_ext, mature_ref, star_ref, species):
		"""Parses the reference files and stores them as class attributes"""
		self.pri_dict = self._parse_fasta(pri_ref)
		self.mature_ext = self._parse_fasta(mature_ext)
		self.star_ext = self._parse_fasta(star_ext)
		self.mature_dict = self._parse_fasta(mature_ref)
		self.star_dict = self._parse_fasta(star_ref)
		self.mature_pos_dict = self._get_sequence_matches(self.mature_dict, self.mature_ext)
		self.star_pos_dict = self._get_sequence_matches(self.star_dict, self.star_ext)


	def _parse_fasta(self, file_path):
		"""Parse FASTA file and return a dictionary"""
		return {key: value for (key, value) in SeqIO.to_dict(SeqIO.parse(open(file_path, 'r'), 'fasta')).items() if key. startswith(self.species)}


	def _get_sequence_matches(self, ref_dict, ext_dict):
		"""Find exact matches between reference and extended sequences"""
		return {key: self.exact_match(str(ref_dict[key].seq), str(ext_dict[key].seq)) for key in ref_dict if key in ext_dict}


	def start_pos(self, RNAME):
		"""Get read start position"""
		gene_name = RNAME.split('_')[0]
		start_5p, start_3p = 999, 999 # set to this value if no annotated star, reconsider
		if gene_name+'_5p' in self.mature_pos_dict.keys() and gene_name+'_3p' in self.mature_pos_dict.keys():  # 5p mature
			start_5p, start_3p = self.mature_pos_dict[gene_name+'_5p']['start'], self.mature_pos_dict[gene_name+'_3p']['start']
		elif gene_name+'_5p' in self.mature_pos_dict.keys() and gene_name+'_3p*' in self.star_pos_dict.keys(): # 3p mature
			start_5p, start_3p = self.mature_pos_dict[gene_name+'_5p']['start'], self.star_pos_dict[gene_name+'_3p*']['start']
		elif gene_name+'_3p' in self.mature_pos_dict.keys() and gene_name+'_5p*' in self.star_pos_dict.keys(): # co-mature
			start_3p, start_5p = self.mature_pos_dict[gene_name+'_3p']['start'], self.star_pos_dict[gene_name+'_5p*']['start']
		return [start_5p, start_3p]


	def end_pos(self, RNAME):
		"""Get read end position"""
		gene_name = RNAME.split('_')[0]
		start_5p, start_3p = 999, 999 # set to this value if no annotated star, reconsider
		if gene_name+'_5p' in self.mature_pos_dict.keys() and gene_name+'_3p' in self.mature_pos_dict.keys():  # 5p mature
			start_5p, start_3p = self.mature_pos_dict[gene_name+'_5p']['end'], self.mature_pos_dict[gene_name+'_3p']['end']
		elif gene_name+'_5p' in self.mature_pos_dict.keys() and gene_name+'_3p*' in self.star_pos_dict.keys(): # 3p mature
			start_5p, start_3p = self.mature_pos_dict[gene_name+'_5p']['end'], self.star_pos_dict[gene_name+'_3p*']['end']
		elif gene_name+'_3p' in self.mature_pos_dict.keys() and gene_name+'_5p*' in self.star_pos_dict.keys(): # co-mature
			start_3p, start_5p = self.mature_pos_dict[gene_name+'_3p']['end'], self.star_pos_dict[gene_name+'_5p*']['end']
		return [start_5p, start_3p]


	def exact_match(self, pattern, text):
		"""Function to identify start and end position of read to reference"""
		pattern = pattern.replace('U', 'T')
		text = text.replace('U', 'T')
		start = text.find(pattern) + 1
		end = start + len(pattern)
		return {'start': start, 'end': end}


	def abs_min(self, x, y, read_side):
		"""Get the absolute minimum value based on the read side"""
		return x if read_side == 'read_5p' else y


	def is_5p_or_3p(self, RNAME):
		"""Check if the read is 5p or 3p"""
		if RNAME.endswith('_5p') or RNAME.endswith('_5p*'):
			return('read_5p')
		if RNAME.endswith('_3p') or RNAME.endswith('_3p*'):
			return('read_3p')


	def get_read_end(self, read_side, canon_5p_start, canon_3p_start, read_start, read_len, canon_5p_end, canon_3p_end):
		"""Calculate the read end position based on the read side"""
		if read_side == 'read_5p':
			return(canon_5p_start + read_start + read_len - canon_5p_end)
		if read_side == 'read_3p':
			return(canon_3p_start + read_start + read_len - canon_3p_end)


	def mature_star(self, RNAME):
		"""Determine if the sequence is a mature miRNA or a star sequence"""
		return 'star' if RNAME.endswith('*') else 'mature'


	def add_counts(self, i, mature_or_star, start_or_end, mature_star, read_start, read_end, counts):
		if read_start == i and mature_star == mature_or_star and start_or_end == 'start':
			return(counts)
		elif read_end == i and mature_star == mature_or_star and start_or_end == 'end':
			return(counts)
		else:
			return(0)


	def get_canon_mature_seq(self, gene_name, read_side, mature_star):
		"""Get canonical mature"""
		if read_side == 'read_5p':
			return str(self.mature_dict[gene_name + '_5p'].seq) if mature_star == 'mature' else str(self.mature_dict[gene_name + '_3p'].seq)
		elif read_side == 'read_3p':
			return str(self.mature_dict[gene_name + '_3p'].seq) if mature_star == 'mature' else str(self.mature_dict[gene_name + '_5p'].seq)


	def get_canon_star_seq(self, gene_name, read_side, mature_star):
		"""Get canonical star sequence"""
		if read_side == 'read_5p':
			return str(self.star_dict[gene_name + '_3p*'].seq) if mature_star == 'mature' else str(self.star_dict[gene_name + '_5p*'].seq)
		elif read_side == 'read_3p':
			return str(self.star_dict[gene_name + '_5p*'].seq) if mature_star == 'mature' else str(self.star_dict[gene_name + '_3p*'].seq)


	def get_canon_type(self,read_side, mature_star, gene_name):
		"""Annotate 5p, 3p or co-mature"""
		type_ = ''
		if read_side == 'read_5p':
			type_ = 'mature_5p' if mature_star == 'mature' else 'mature_3p'
		elif read_side == 'read_3p':
			type_ = 'mature_3p' if mature_star == 'mature' else 'mature_5p'
		if gene_name + '_5p' in self.mature_dict and gene_name + '_3p' in self.mature_dict:
			type_ = 'co_mature'
		return type_


	def add_zero_rows(self, tissue):
		non_counted_genes = [gene_name[:-4] for gene_name in self.pri_dict.keys() if gene_name[:-4] not in self.df.Gene_name.values]
		if len(non_counted_genes) > 0:
			for i in non_counted_genes:
				if i+'_5p' in self.mature_dict.keys():
					mature_seq = str(self.mature_dict[i+'_5p'].seq)
					mature_type = 'mature_5p'
				if i+'_3p' in self.mature_dict.keys():
					mature_seq = str(self.mature_dict[i+'_3p'].seq)
					mature_type = 'mature_3p'
				if i+'_5p' in self.mature_dict.keys() and i+'_3p' in self.mature_dict.keys():
					mature_type = 'co_mature'
				if i+'_3p*' in self.star_dict.keys():
					star_seq = str(self.star_dict[i+'_3p*'].seq)
				elif i+'_5p*' in self.star_dict.keys():
					star_seq = str(self.star_dict[i+'_5p*'].seq)
				else:
					star_seq = ''
				self.df = self.df.append({
						'Gene_name' : i, 'Mature_Seq' : mature_seq, 'Star_Seq' : star_seq, 
						'Mature_type' : mature_type, 'Tissue' : tissue, 'Species' : species_names[self.species],
						'Mature_Start_-5': 0 , 'Mature_Start_-4': 0 , 'Mature_Start_-3': 0 , 'Mature_Start_-2': 0 , 
						'Mature_Start_-1': 0 , 'Mature_Start_0': 0 , 'Mature_Start_1': 0 , 'Mature_Start_2': 0 , 
						'Mature_Start_3': 0 , 'Mature_Start_4': 0 , 'Mature_Start_5': 0 , 'Mature_end_-5': 0 , 
						'Mature_end_-4': 0 , 'Mature_end_-3': 0 , 'Mature_end_-2': 0 , 'Mature_end_-1': 0 , 
						'Mature_end_0': 0 , 'Mature_end_1': 0 , 'Mature_end_2': 0 , 'Mature_end_3': 0 , 
						'Mature_end_4': 0 , 'Mature_end_5': 0 , 'Star_Start_-5': 0 , 'Star_Start_-4': 0 , 
						'Star_Start_-3': 0 , 'Star_Start_-2': 0 , 'Star_Start_-1': 0 , 'Star_Start_0': 0 , 
						'Star_Start_1': 0 , 'Star_Start_2': 0 , 'Star_Start_3': 0 , 'Star_Start_4': 0 , 'Star_Start_5': 0 , 
						'Star_End_-5': 0 , 'Star_End_-4': 0 , 'Star_End_-3': 0 , 'Star_End_-2': 0 , 'Star_End_-1': 0 ,
						'Star_End_0': 0 , 'Star_End_1': 0 , 'Star_End_2': 0 , 'Star_End_3': 0 , 'Star_End_4': 0 ,
						'Star_End_5': 0 
						}, ignore_index=True)


	def count_total_reads(self, mature_star):
		colnames = [mature_star+'_Start_'+str(i) for i in range(-5, 5 + 1)]
		return(self.df[colnames].sum(axis=1))


	def sam2df(self, file_path):
		"""
		Convert a samfile into a pandas dataframe
		"""
		with open(file_path) as f:
			"""Skip the header lines"""
			skiprows = 0
			for line in f:
				if line.startswith('@'):
					skiprows += 1
				else:
					break
			f.close()
		"""Read sam file as a tsv file, skip the header lines"""
		self.df = pd.read_csv(file_path, sep='\t', skiprows=skiprows, 
			names = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 
					   'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'OPT1',
					   'OPT2', 'OPT3'], warn_bad_lines=False, error_bad_lines=False, low_memory=False)
		self.df = self.df[self.df.RNAME != '*']
		self.df['read_len']       = [int(x.split('M')[0]) for x in self.df.CIGAR]
		self.df['canon_5p_start'] = [self.start_pos(RNAME)[0] for RNAME in self.df.RNAME]
		self.df['canon_3p_start'] = [self.start_pos(RNAME)[1] for RNAME in self.df.RNAME]
		self.df['canon_5p_end']   = [self.end_pos(RNAME)[0] for RNAME in self.df.RNAME]
		self.df['canon_3p_end']   = [self.end_pos(RNAME)[1] for RNAME in self.df.RNAME]
		self.df['diff_5p_start']  = self.df.POS - self.df.canon_5p_start
		self.df['diff_3p_start']  = self.df.POS - self.df.canon_3p_start
		self.df['read_side']      = list(map(self.is_5p_or_3p, self.df.RNAME))
		self.df['read_start']     = list(map(self.abs_min, self.df.diff_5p_start, self.df.diff_3p_start, self.df.read_side))
		self.df['read_end']       = list(map(self.get_read_end, self.df.read_side, self.df.canon_5p_start, self.df.canon_3p_start, self.df.read_start, self.df.read_len, self.df.canon_5p_end, self.df.canon_3p_end))
		self.df['mature_star']    = list(map(self.mature_star, self.df.RNAME))
		self.df['counts']         = [int(x.split('-')[1]) for x in self.df.QNAME]
		self.df = self.df[(self.df.read_start <= 5) & (self.df.read_start >= -5)]
		for i in range(-5, 5+1):
			self.df['Mature_Start_' + str(i)] = list(map(self.add_counts, [i]*len(self.df), ['mature']*len(self.df), ['start']*len(self.df), self.df.mature_star, self.df.read_start, self.df.read_end, self.df.counts))
		for i in range(-5, 5+1):	
			self.df['Mature_end_'   + str(i)] = list(map(self.add_counts, [i]*len(self.df), ['mature']*len(self.df), ['end']*len(self.df), self.df.mature_star, self.df.read_start, self.df.read_end, self.df.counts))
		for i in range(-5, 5+1):
			self.df['Star_Start_'   + str(i)] = list(map(self.add_counts, [i]*len(self.df), ['star']*len(self.df),   ['start']*len(self.df), self.df.mature_star, self.df.read_start, self.df.read_end, self.df.counts)) 
		for i in range(-5, 5+1):
			self.df['Star_End_'     + str(i)] = list(map(self.add_counts, [i]*len(self.df), ['star']*len(self.df),   ['end']*len(self.df),   self.df.mature_star, self.df.read_start, self.df.read_end, self.df.counts))
		self.df['Gene_name'] = [RNAME.split('_')[0] for RNAME in self.df.RNAME]
		self.df['Mature_Seq']  = list(map(self.get_canon_mature_seq, self.df.Gene_name, self.df.read_side, self.df.mature_star))
		self.df['Star_Seq']    = list(map(self.get_canon_star_seq, self.df.Gene_name, self.df.read_side, self.df.mature_star))
		self.df['Tissue']    = file_path.split('/')[-1]
		self.df['Species']   = species_names[self.species]
		self.df['Mature_type'] = list(map(self.get_canon_type, self.df.read_side, self.df.mature_star, self.df.Gene_name))
		self.df = self.df[['Gene_name', 'Mature_Seq', 'Star_Seq', 'Mature_type', 'Mature_Start_-5', 'Mature_Start_-4', 'Mature_Start_-3', 'Mature_Start_-2', 'Mature_Start_-1',
		'Mature_Start_0', 'Mature_Start_1', 'Mature_Start_2', 'Mature_Start_3', 'Mature_Start_4', 'Mature_Start_5', 'Mature_end_-5', 'Mature_end_-4',
		'Mature_end_-3', 'Mature_end_-2', 'Mature_end_-1', 'Mature_end_0', 'Mature_end_1', 'Mature_end_2', 'Mature_end_3', 'Mature_end_4',
		'Mature_end_5', 'Star_Start_-5', 'Star_Start_-4', 'Star_Start_-3', 'Star_Start_-2', 'Star_Start_-1', 'Star_Start_0', 'Star_Start_1',
		'Star_Start_2', 'Star_Start_3', 'Star_Start_4', 'Star_Start_5', 'Star_End_-5', 'Star_End_-4', 'Star_End_-3', 'Star_End_-2',
		'Star_End_-1', 'Star_End_0', 'Star_End_1', 'Star_End_2', 'Star_End_3', 'Star_End_4', 'Star_End_5', 'Tissue', 'Species']]
		pd.set_option('display.max_columns', None)
		print(self.df[self.df.Gene_name == 'Hsa-Mir-154-P28b'])
		self.df = self.df.groupby(by=['Gene_name', 'Mature_Seq', 'Mature_type', 'Tissue', 'Species']).sum()
		self.df = self.df.reset_index()
		print(self.df[self.df.Gene_name == 'Hsa-Mir-154-P28b'])

		self.add_zero_rows(file_path.split('/')[-1])
		self.df['Mature_total'] = self.count_total_reads('Mature')
		self.df['Star_total']   = self.count_total_reads('Star')


# Loop over sams 
# ### Should instead be snakemake compatible script running individual samfiles, followed by aggregating rule.

samdir = 'mirgenedb_merged_extended_sam_filtered_no_pos2_18_mismatch/'
outdir = 'aligned_positions_extended_mature/'
skipfiles = ['Muscle_Biceps_female_homogenized.sam', 'female_soma_SRR5602511.sam']

for directory in os.listdir(samdir):

    if directory != 'hsa':
        continue

    if not os.path.isdir(samdir + str(directory)):
        continue
    species_id = str(directory).title()
    species = species_names[species_id]

    print(species_id)
    
    if os.path.isdir(outdir+str(directory)):
        for file in os.listdir(outdir+directory):
            if file.endswith('.csv'):
                os.remove(outdir+directory+'/'+file)
    
    species_obj = Correlomics(species_id)
    print('nr_files: ',len(glob.glob(samdir + str(directory) + '/*.sam')))
    if len(glob.glob(samdir + str(directory) + '/*.sam')) == 0:
        species_obj.df = pd.DataFrame(columns=['Gene_name', 'Mature_Seq', 'Star_Seq', 'Mature_type', 'Tissue', 'Species', 'Mature_Start_-5', 'Mature_Start_-4', 
		'Mature_Start_-3', 'Mature_Start_-2', 'Mature_Start_-1',
            'Mature_Start_0', 'Mature_Start_1', 'Mature_Start_2', 'Mature_Start_3', 'Mature_Start_4', 'Mature_Start_5', 'Mature_end_-5', 'Mature_end_-4',
            'Mature_end_-3', 'Mature_end_-2', 'Mature_end_-1', 'Mature_end_0', 'Mature_end_1', 'Mature_end_2', 'Mature_end_3', 'Mature_end_4',
            'Mature_end_5', 'Star_Start_-5', 'Star_Start_-4', 'Star_Start_-3', 'Star_Start_-2', 'Star_Start_-1', 'Star_Start_0', 'Star_Start_1',
            'Star_Start_2', 'Star_Start_3', 'Star_Start_4', 'Star_Start_5', 'Star_End_-5', 'Star_End_-4', 'Star_End_-3', 'Star_End_-2',
            'Star_End_-1', 'Star_End_0', 'Star_End_1', 'Star_End_2', 'Star_End_3', 'Star_End_4', 'Star_End_5'])
        species_obj.df.Tissue = 'total_animal.sam'
        species_obj.add_zero_rows('total_animal.sam')
        if not os.path.isdir(outdir+str(directory)+'/'):
            os.mkdir(outdir+str(directory)+'/')
        species_obj.df['Mature_total'] = species_obj.count_total_reads('Mature')
        species_obj.df['Star_total'] = species_obj.count_total_reads('Star')
        species_obj.df.to_csv(outdir+str(directory)+'/'+'aligned_positions_whole_animal.csv')
    for filename in os.listdir(samdir + str(directory)):
        if filename.startswith('.'):
            continue
        if not filename.endswith('.sam'):
            continue
        if filename in skipfiles:
            continue
        print('    ', filename)
        species_obj.sam2df(samdir+str(directory+'/'+str(filename)))
        if not os.path.isdir(outdir+str(directory)+'/'):
            os.mkdir(outdir+str(directory)+'/')
        species_obj.df.to_csv(outdir+str(directory)+'/'+'aligned_positions_'+filename.split('.sam')[0]+'.csv')

path = '/home/jcdenton/projects/mirgenedb_editing_events/aligned_positions_extended_mature/'
li = []
for directory in os.listdir(path):
    all_files = glob.glob(path+directory+'/*.csv')
    for filename in all_files:
        df = pd.read_csv(filename)
        li.append(df)
frame = pd.concat(li, axis=0, ignore_index=True)
frame.to_csv(path+'merged_aligned_positions_extended_mature.csv')

