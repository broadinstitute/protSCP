import argparse
from .preprocess_mutatex import preprocess_mutatex_data
from .preprocess_structure import preprocess_structure
from .merge_datasets import merge_protein_datasets


""" 
Function: prepare_protein_dataset
Description: Preprocess raw mutatex results and Protein Structure data.
Args:
    [required] raw_mutatex_filename: path to raw mutatex results
    [required for AlphaFoldDB Structure] uniprot_id: UniProt ID of the protein
    [required for mmCIF Structure] cif_path: path to mmCIF file
Returns:
    cluster_df: pandas dataframe with mutatex and structure data for clustering
"""
def prepare_protein_dataset(raw_mutatex_filename, uniprot_id=None, cif_path=None):
    mutatex_df = preprocess_mutatex_data(raw_mutatex_filename)
    structure_df = preprocess_structure(uniprot_id, cif_path)
    merged = merge_protein_datasets([mutatex_df, structure_df])
    return merged

def main():
    parser = argparse.ArgumentParser(description='Prepare protein dataset')
    parser.add_argument('raw_mutatex_filename', type=str, help='Raw mutatex filename')
    parser.add_argument('uniprot_id', type=str, help='UniProt ID of the protein')
    parser.add_argument('output_file', type=str, help='Output file')
    args = parser.parse_args()

    df = prepare_protein_dataset(args.raw_mutatex_filename, args.uniprot_id)
    df.to_csv(args.output_file, index=False)

if __name__ == "__main__":
    main()
