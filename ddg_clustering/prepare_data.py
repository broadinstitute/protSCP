import argparse
from .preprocess_mutatex import preprocess_mutatex_data
from .preprocess_structure import preprocess_structure
from .merge_datasets import merge_protein_datasets

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
