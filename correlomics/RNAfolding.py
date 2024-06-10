import os
from Bio import SeqIO
import pickle

# Define file paths
BASE_PATH = '/home/jcdenton/projects/mirgenedb_database/'
STEM_FILE_PATH = os.path.join(BASE_PATH, 'merged_all-pri-30-30.fas')
MATURE_FILE_PATH = os.path.join(BASE_PATH, 'mature.fas')
STAR_FILE_PATH = os.path.join(BASE_PATH, 'star.fas')

# Open files
with open(STEM_FILE_PATH, 'r') as stem_file, \
        open(MATURE_FILE_PATH, 'r') as mature_file, \
        open(STAR_FILE_PATH, 'r') as star_file:
    # Parse fasta files
    stem = SeqIO.to_dict(SeqIO.parse(stem_file, "fasta"))
    mature = SeqIO.to_dict(SeqIO.parse(mature_file, "fasta"))
    star = SeqIO.to_dict(SeqIO.parse(star_file, "fasta"))


# mfold environment variables
MFOLD_PATH = '/home/jcdenton/projects/mirgenedb_editing_events__deprecated__/src/RNAstructure/'
os.environ['PATH'] += os.pathsep + MFOLD_PATH + 'exe'
os.environ['DATAPATH'] = MFOLD_PATH + 'data_tables/'


class RNAfold():
    def __init__(self, mirgenedb_id):
        """
        Initialize the Fold object.


        Args:
        - mirgenedb_id (str): The identifier of the RNA sequence.
        """
        self.mirgenedb_id = mirgenedb_id
        self.stem, self.mature, self.star, self._5p_arm, self._3p_arm = self.fetch_sequences()

    def is_5p_or_3p_mature(self):
        """
        Determine if the RNA is 5' or 3' mature.

        Returns:
        - str: Either "co_mature", "5p_mature", or "3p_mature".
        """
        mature_id = [id for id in mature]
        if self.mirgenedb_id + "_5p" in mature_id and self.mirgenedb_id + "_3p" in mature_id:
            return "co_mature"
        elif self.mirgenedb_id + "_5p" in mature_id:
            return "5p_mature"
        else:
            return "3p_mature"

    def fetch_sequences(self):
        """
        Fetch RNA sequences from files based on maturity.

        Returns:
        - tuple: A tuple containing the stem, mature, star, 5' arm, and 3' arm
        sequences of the RNA.
        """
        mature_arm = self.is_5p_or_3p_mature()
        if mature_arm == "5p_mature":
            return (str(stem[self.mirgenedb_id + "_pri"].seq), 
                    str(mature[self.mirgenedb_id + "_5p"].seq), 
                    str(star[self.mirgenedb_id + "_3p*"].seq),
                    str(mature[self.mirgenedb_id + "_5p"].seq),
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

    def exact_match(self, pattern, text):
        """
        Identify start and end position of a pattern in a text.

        Args:
        - pattern (str): The pattern to search for.
        - text (str): The text to search within.

        Returns:
        - dict: A dictionary containing the start and end positions of the
        pattern.
        """
        pattern = pattern.replace('U', 'T')
        text = text.replace('U', 'T')
        start = text.find(pattern) + 1
        end = start + len(pattern)
        return {'start': start, 'end': end}

    def write_constraints(self, constrain_stem, constrain_loop, no_pairing, sequence_length, midpoint):
        if no_pairing:
            prohibit_pairs = [
                f"P {1} {0} {midpoint}",
                f"P {midpoint + 1} {0} {sequence_length}"
            ]
            with open('constraints.txt', 'w') as constraints:
                for pair in prohibit_pairs:
                    constraints.write(f"{pair}\n")
        else:
            with open('constraints.txt', 'w') as constraints:
                if constrain_stem:
                    constraints.write(f'P {1} {0} {17}\n')  # Define 5p ds
                    constraints.write(f'P {len(self.stem) - 17 + 1} {0} {17}\n')  # Define 3p ds
                if constrain_loop:
                    constraints.write(f'P {self.end_5p_arm + 1} {0} {self.start_3p_arm - self.end_5p_arm - 2}\n')

    def fold_hairpin(self, constrain_stem=False, constrain_loop=False, no_pairing=False):
        """
        Fold the RNA hairpin structure and extract relevant information.
        """
        self.end_5p_arm = self.exact_match(self._5p_arm, self.stem)['end']
        self.start_3p_arm = self.exact_match(self._3p_arm, self.stem)['start']
        sequence_length = len(self.stem)
        midpoint = sequence_length // 2

        os.chdir('tmp')

        with open('workingfile.txt', 'w') as workingfile:
            workingfile.write(self.stem + '\n')

        self.write_constraints(constrain_stem, constrain_loop, no_pairing, sequence_length, midpoint)

        os.system('mfold SEQ=workingfile.txt AUX=constraints.txt')
        os.system('ct2dot workingfile.txt.ct ALL workingfile.txt.dot')

        with open('workingfile.txt.dot', 'r') as file:
            self.dG = file.readline()
            self.seq = file.readline()
            self.dot = file.readline()

        

        os.chdir('..')

    def find_motif(self):
        # Determine positions of 5' and 3' arms
        self.end_5p_arm = self.exact_match(self._5p_arm, self.stem)['end']
        self.start_3p_arm = self.exact_match(self._3p_arm, self.stem)['start']

        # Calculate positions relative to arms
        pos_5p = [i + self.end_5p_arm for i in [-4, -3, -2, -1]]
        pos_3p = [i + self.start_3p_arm for i in [1, 0, -1, -2]]
        
        # Extract sequences and dot-bracket notations
        dicer_cleave = [
            [self.seq[i] for i in pos_5p],
            [self.dot[i] for i in pos_5p],
            [self.dot[i] for i in pos_3p],
            [self.seq[i] for i in pos_3p]
                        ]

        # Initialize motif list
        self.GYM = [0, 0, 0]

        # Check for specific motifs related to dicer cleavage
        if dicer_cleave[3][-2] == 'G':
            self.GYM[0] = 1
        if dicer_cleave[3][-3] in ('C', 'T', 'U'):
            self.GYM[1] = 1
        if dicer_cleave[2][-4] == '.':
            self.GYM[2] = 1

    def count_bulges(self, fold):
        return fold.count('.')

    def count_basepair(self, fold):
        return sum(1 for i in fold if i in ('(', ')'))

    def check_misfold(self):
        
        # Determine positions of 5' and 3' arms
        self.end_5p_arm = self.exact_match(self._5p_arm, self.stem)['end']
        self.start_3p_arm = self.exact_match(self._3p_arm, self.stem)['start']

        # Check if brackets are wrong direction upstream or downstream 
        for i, dot in enumerate(self.dot):
            if (dot == ')' and i < self.end_5p_arm) or (dot == '(' and i > self.start_3p_arm):
                return True
        return False

    def find_bulge_length(self):
        # Determine positions of 5' and 3' arms
        pos_5p_arm = self.exact_match(self._5p_arm, self.stem)
        pos_3p_arm = self.exact_match(self._3p_arm, self.stem)
    
        # Extract dot-bracket notation for arms
        bulges_5p = self.dot[pos_5p_arm['start']:pos_5p_arm['end']]
        bulges_3p = self.dot[pos_3p_arm['start']:pos_3p_arm['end']-2] # Get correct counting

        # Calculate bulge length and characteristics
        self.bulge_len_dict = {
                '5p_length': len(self._5p_arm),
                '5p_bulges': self.count_bulges(bulges_5p),
                '5p_bulges_len': float(self.count_bulges(bulges_5p)) / (len(self._5p_arm)),
                '5p_basepair': self.count_basepair(bulges_5p),
                '5p_basepair_len': float(self.count_basepair(bulges_5p)) / (len(self._5p_arm)),
                '3p_length': len(self._3p_arm),
                '3p_bulges': self.count_bulges(bulges_3p),
                '3p_bulges_len': float(self.count_bulges(bulges_3p)) / len(self._3p_arm),
                '3p_basepair': self.count_basepair(bulges_3p),
                '3p_basepair_len': float(self.count_basepair(bulges_3p)) / (len(self._3p_arm)) 
                }

    def fetch_hairpin(self, path):
        with open(path, 'rb') as file:
            self.hairpin = pickle.load(file)
            self.seq = self.hairpin[1]
            self.dot = self.hairpin[2]
            self.dG = self.hairpin[3]

