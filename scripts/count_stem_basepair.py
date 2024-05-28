import os
import sys
import pandas as pd

# Add the parent directory of your current script to the Python path
current_dir = os.path.dirname(__file__)
parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
sys.path.append(parent_dir)

# Now you can import your module
from correlomics.foldings import Fold

hairpin_dict = '/home/jcdenton/projects/correlomics/test_folds/'
bulge_len_dict = {}
misfolded_ids = []

# Open file for writing misfolds
with open('misfolds.txt', 'w') as misfold_file:
    for file in os.listdir(hairpin_dict):
        mirgenedb_id = file[:-4]
        fold = Fold(mirgenedb_id)
        fold.fetch_hairpin(hairpin_dict + file)

        # Check for misfold and write details to file
        if fold.check_misfold():
            misfolded_ids.append(mirgenedb_id)
            misfold_file.write(f"{mirgenedb_id}\n")
            misfold_file.write(f"{fold.hairpin['seq']}\n")
            misfold_file.write(f"{fold.hairpin['dot']}\n")
            continue

        # Find bulge length and store in dictionary
        fold.find_bulge_length()
        bulge_len_dict[mirgenedb_id] = fold.bulge_len_dict

# Convert dictionary to DataFrame and save to CSV
df = pd.DataFrame.from_dict(bulge_len_dict, orient='index')
column_names = ['5p_length', '5p_bulges', '5p_bulges_len', '5p_basepair',
                '5p_basepair_len', '3p_length', '3p_bulges', '3p_bulges_len',
                '3p_basepair', '3p_basepair_len']
df.columns = column_names
df.to_csv('bulge_len_ratio.csv')
