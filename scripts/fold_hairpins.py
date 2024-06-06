import os
import sys
import yaml
import pickle
import argparse
import pandas as pd
from Bio import SeqIO

current_dir = os.path.dirname(__file__)
parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
sys.path.append(parent_dir)
from correlomics.foldings import Fold


def load_config(config_file):
    """Load configuration from a YAML file."""
    with open(config_file, 'r') as file:
        config = yaml.safe_load(file)
        return config


def load_sequences(file_paths):
    """Load sequences from given file paths."""
    sequences = {}
    for key, file_path in file_paths.items():
        with open(file_path, 'r') as file:
            sequences[key] = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
    return sequences


def save_hairpin(hairpin, output_dir, suffix):
    os.makedirs(output_dir, exist_ok=True)  # Ensure the output directory exists
    file_path = os.path.join(output_dir, f"{hairpin.mirgenedb_id}.pkl")
    with open(file_path, 'wb') as f:
        pickle.dump([hairpin.mirgenedb_id+'\n', hairpin.seq, hairpin.dot, hairpin.dG], f)


def process_hairpins(sequences, output_dirs, filters):
    """Process hairpins and save them to respective directories."""
    misfold_dict = {}
    exclude_list = filters.get('exclude', [])
    include_only_list = filters.get('include_only', [])

    for id in sequences['pri_mirna'].keys():
        if any(excl in id for excl in exclude_list):
            continue
        if include_only_list and not any(inc in id for inc in include_only_list):
            continue

        misfold_dict[id] = [0, 0, 0, 0]
        constraints = [
            (False, False, 'none'),
            (True, False, 'stem'),
            (True, True, 'both'),
            (False, True, 'loop')
        ]

        for i, (constrain_stem, constrain_loop, suffix) in enumerate(constraints):
            hairpin = Fold(id[:-4])
            hairpin.fold_hairpin(constrain_stem=constrain_stem, constrain_loop=constrain_loop)
            save_hairpin(hairpin, output_dirs[suffix], suffix)
            if hairpin.check_misfold():
                misfold_dict[id][i] = 1

    return misfold_dict


def main():
    # Define command line arguments
    parser = argparse.ArgumentParser(description='Fold pri-miRNA hairpin structures')
    parser.add_argument('config_file', type=str, help='Path to YAML config')
    args = parser.parse_args()

    # Load configuration from YAML file
    config = load_config(args.config_file)
    base_path = config['file_paths']['base_path']
    file_paths = {
        key: os.path.join(base_path, relative_path)
        for key, relative_path in config['file_paths']['fasta_files'].items()
    }

    sequences = load_sequences(file_paths)

    output_dirs = {
        key: relative_path
        for key, relative_path in config['out_dir'].items()
    }

    filters = config.get('filters', {})
    misfold_dict = process_hairpins(sequences, output_dirs, filters)

    df = pd.DataFrame.from_dict(misfold_dict, orient='index', columns=['none', 'stem', 'loop', 'both'])
    print(df)


if __name__ == "__main__":
    main()
