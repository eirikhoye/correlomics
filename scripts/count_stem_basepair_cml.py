import os
import sys
import pandas as pd
import argparse

# Define command line arguments
parser = argparse.ArgumentParser(description='Count basepairs in hairpin stem')
parser.add_argument('hairpin_dir', type=str, help='Directory containing hairpin files')
parser.add_argument('--output_csv', type=str, default='bulge_len_ratio.csv', help='Output CSV file for bulge lengths')

args = parser.parse_args()

# Add the parent directory of your current script to the Python path
current_dir = os.path.dirname(__file__)
parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
sys.path.append(parent_dir)

# Now you can import your module
from correlomics.RNAfolding import Fold

bulge_len_dict = {}
misfolded_ids = []

# Open file for writing misfolds
for file in os.listdir(args.hairpin_dir):
    mirgenedb_id = file[:-4]
    fold = Fold(mirgenedb_id)
    fold.fetch_hairpin(os.path.join(args.hairpin_dir, file))
    
    # Find bulge length and store in dictionary
    fold.find_bulge_length()
    bulge_len_dict[mirgenedb_id] = fold.bulge_len_dict

# Convert dictionary to DataFrame and save to CSV
df = pd.DataFrame.from_dict(bulge_len_dict, orient='index')
column_names = ['5p_length', '5p_bulges', '5p_bulges_len', '5p_basepair',
                '5p_basepair_len', '3p_length', '3p_bulges', '3p_bulges_len',
                '3p_basepair', '3p_basepair_len']

df.columns = column_names
df.to_csv(args.output_csv)

