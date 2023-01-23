"""
A potential refactor of the Correlomics class could be to extract the file parsing logic into a separate method, and to extract the logic for finding the start and end positions of a given gene into separate methods as well. Additionally, it would be good to add some docstrings to explain what the class and methods are doing. Also, it would be good to move all the import statements inside the class and remove the unnecessary imports.
Here is an example refactor:

"""

import json
from Bio import SeqIO

class Correlomics():
    def __init__(self, species):
        """Initialize the class with the species name, and parse the reference files"""
        self.species = species
        with open('species_dict.txt') as f:
            data = f.read()
        self.species_names = json.loads(data)
        self.parse_reference_files('data/fasta/prifile.fasta', 'data/fasta/mature_extended.fas', 'data/fasta/star_extended_fake_451.fas', 'data/fasta/mature.fasta', 'data/fasta/starfile_fake_451.fasta', self.species)

    def parse_reference_files(self, pri_ref, mature_ext, star_ext, mature_ref, star_ref, species):
        """Parses the reference files and stores them in class attributes"""
        self.pri_dict        = {key:value for (key,value) in SeqIO.to_dict(SeqIO.parse(open(pri_ref,    'r'), 'fasta')).items() if key.startswith(self.species)}
        self.mature_ext      = {key:value for (key,value) in SeqIO.to_dict(SeqIO.parse(open(mature_ext, 'r'), 'fasta')).items() if key.startswith(self.species)}
        self.star_ext        = {key:value for (key,value) in SeqIO.to_dict(SeqIO.parse(open(star_ext,   'r'), 'fasta')).items() if key.startswith(self.species)}
        self.mature_dict     = {key:value for (key,value) in SeqIO.to_dict(SeqIO.parse(open(mature_ref, 'r'), 'fasta')).items() if key.startswith(self.species)}
        self.star_dict       = {key:value for (key,value) in SeqIO.to_dict(SeqIO.parse(open(star_ref,   'r'), 'fasta')).items() if key.startswith(self.species)}
        self.mature_pos_dict = {key : self.exact_match(str(self.mature_dict[key].seq), str(self.mature_ext[key].seq)) for key, value in self.mature_dict.items() if key in self.mature_ext}
        self.star_pos_dict   = {key : self.exact_match(str(self.star_dict[key].seq),   str(self.star_ext[key].seq)) for key, value in self.star_dict.items()
