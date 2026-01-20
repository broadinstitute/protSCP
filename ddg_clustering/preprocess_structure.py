import argparse
import requests
import pandas as pd
from Bio.PDB import MMCIFParser
from io import StringIO
# get 3 to 1 map from Bio.PDB
from Bio.Data.PDBData import protein_letters_3to1


def download_structure(uniprot_id):
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v6.cif"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        raise RuntimeError("Failed to download structure. Please check the UniProt ID and try again.")

def parse_structure(cif_data):
    parser = MMCIFParser()
    file_like = StringIO(cif_data)
    structure = parser.get_structure('protein', file_like)
    model = structure[0]

    data = []
    for chain in model:
        for residue in chain:
            if 'CA' in residue:
                ca = residue['CA']
                # get plddt from B-factor
                plddt = ca.get_bfactor()
                aa3 = residue.resname
                # convert 3-letter residue name to 1-letter
                aa1 = protein_letters_3to1.get(aa3, 'X')

                data.append([
                    aa1, chain.id, residue.id[1],
                    ca.coord[0], ca.coord[1], ca.coord[2], plddt
                ])

    df = pd.DataFrame(data, columns=['aa', 'chain', 'residue_id', 'xca', 'yca', 'zca', 'plddt'])
    return df

def preprocess_structure(uniprot_id=None, cif_path=None):
    if cif_path is not None:
        with open(cif_path, 'r') as file:
            cif_data = file.read()
            df = parse_structure(cif_data)
    elif uniprot_id is not None:
        cif_data = download_structure(uniprot_id)
        df = parse_structure(cif_data)
    else:
        raise ValueError("Either cif_path or uniprot_id must be provided")
      
    return df

def main():
    parser = argparse.ArgumentParser(description='Download and parse AlphaFold protein structure.')
    parser.add_argument('uniprot_id', type=str, help='UniProt ID of the protein')
    parser.add_argument('output_file', type=str, help='Output CSV file')
    args = parser.parse_args()

    df = preprocess_structure(args.uniprot_id)
    df.to_csv(args.output_file, index=False)

if __name__ == "__main__":
    main()
