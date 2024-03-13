import pandas as pd
import argparse


def load_data(input_file):
    return pd.read_csv(input_file, sep='\t')


def split_name(mismatch):
    if pd.isna(mismatch):
        mismatch = 'None_None'
    return str(mismatch).split('_')


def split_edit(mismatch):
    edit_from, edit_to = str(mismatch).split('to')
    return edit_from, edit_to.split('_')[0]


def process_data(df, columns):
    group_cols = ['miRNA'] + (['tissue'] if columns == 'tissue' else [])

    # Calculate total reads per miRNA
    df['total'] = df.groupby(group_cols)['nr_reads'].transform('sum')

    # Convert from wide to long, and drop unnecessary columns
    df = df.melt(id_vars=group_cols + ['ID', 'nr_reads', 'arm',
                        'mature_star', 'read_start', 'read_seq',
                        'total'],
                    value_vars=['from_to_1', 'from_to_2', 'from_to_3'],
                    value_name='mismatch').drop('variable', axis=1)

    # Calculate total number of mismatches
    total_mismatch = df.groupby(
        group_cols + ['mismatch'])['nr_reads'].transform('sum')\
        .rename('total_mismatch').reset_index()

    df['total_mismatch'] = total_mismatch['total_mismatch']

    # Select relevant columns, drop duplicates and sort the DataFrame
    df = df[group_cols + ['arm', 'mature_star', 'total',
                          'mismatch', 'total_mismatch']] \
        .drop_duplicates() \
        .dropna() \
        .sort_values(by=['total', 'total_mismatch'], ascending=False)

    # Reset index
    df.reset_index(inplace=True, drop=True)

    # Split the mismatch name and extract relevant information
    df['mismatch_pos'] = df['mismatch'].apply(lambda x: split_name(x)[1])
    df['edit_from'] = df['mismatch'].apply(lambda x: split_edit(x)[0])
    df['edit_to'] = df['mismatch'].apply(lambda x: split_edit(x)[1])
    df['mismatch'] = df['mismatch'].apply(lambda x: split_name(x)[0])

    # Rearrange and return the DataFrame
    return df[group_cols + ['arm', 'mature_star', 'total',
               'mismatch', 'edit_from', 'edit_to', 'mismatch_pos',
               'total_mismatch']]


def write_output(df, out_file):
    df.to_csv(out_file, sep='\t', index=False)


def main():
    parser = argparse.ArgumentParser(
        description="""
        Process mismatch read files into one editing event per miRNA position.
        """)
    parser.add_argument("input_file", help="Path to the input file")
    parser.add_argument("out_file", help="Path to the output file")
    parser.add_argument("species_tissue", help="Collapse to species or tissue")
    args = parser.parse_args()

    # Check if all command-line arguments are present
    if None in (args.input_file, args.out_file, args.species_tissue):
        parser.error("Missing one or more required arguments")
    if args.species_tissue not in ['species', 'tissue']:
        parser.error("Specify species or tissue.")

    # Run script
    df = load_data(args.input_file)
    processed_df = process_data(df, columns=args.species_tissue)
    write_output(processed_df, args.out_file)


if __name__ == "__main__":
    main()
