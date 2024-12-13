import pandas as pd
import argparse

def merge_protein_datasets(input_dfs):
    """
    Merge protein datasets
    """
    merge_cols = ['residue_id', 'aa', 'chain']

    
    # Merge datasets
    merged = input_dfs[0]
    for dataset in input_dfs[1:]:
        merged = pd.merge(merged, dataset, on=merge_cols, how='outer')

    return merged
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge protein datasets')
    parser.add_argument('input_files', nargs='+', help='Input files')
    parser.add_argument('output_file', help='Output file')
    args = parser.parse_args()

    input_dfs = [pd.read_csv(f) for f in args.input_files]
    merged = merge_protein_datasets(input_dfs)
    merged.to_csv(args.output_file, index=False)
