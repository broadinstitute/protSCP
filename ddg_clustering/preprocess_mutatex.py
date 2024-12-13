import os
import pandas as pd


def split_filename(ddg_file):
    return ddg_file[0], ddg_file[1], int(ddg_file[2:])

def ddg_row_from_file(ddg_results_path):
    def ddg_row_from_file_curried(ddg_file) -> pd.Series:
        residue_type, chain_id, residue_index = split_filename(ddg_file)

        ddg_output_csv = pd.read_csv(ddg_results_path + os.sep + ddg_file, sep=" |\t", engine="python")
        ddg_output_csv.columns = ddg_output_csv.columns[1:].append(pd.Index([0])) # shift column names over 1

        ddg_row_data = ddg_output_csv['avg']
        ddg_row_data.name = (residue_type, chain_id, residue_index)

        return ddg_row_data

    return ddg_row_from_file_curried


def preprocess_mutatex_data(raw_results_path):
    """
    Convert mutatex files from one file per residue to a single file with one row per residue.
    For each residue, there will be one column per mutagenisis run with the avg DDG value from the run.
    """
    ddg_results_path = f"{raw_results_path}/results/mutation_ddgs/final_averages"
    mutation_list_path = f"{raw_results_path}/mutation_list.txt"

    ddg_files = os.listdir(ddg_results_path)

    ddg_rows = list(map(ddg_row_from_file(ddg_results_path), ddg_files)) # why are we sorting .sort(key=residue_id, chain_id)

    mutation_list = pd.read_csv(mutation_list_path, header=None)
    ddg_df = pd.concat(ddg_rows, axis=1).T
    ddg_df.columns = mutation_list[0]

    ddg_df.index.names = ["aa", "chain", "residue_id"]
    
    return ddg_df.sort_values(by=["chain","residue_id"])

def preprocess_mutatex(raw_results_path, output_path):
    mutatex_df = preprocess_mutatex_data(raw_results_path)

    mutatex_df.to_csv(output_path)

def run_multi(inputs_dir, outputs_dir):
    import os
    inputs = os.listdir(inputs_dir)

    # create output dir if it doesn't exist
    if not os.path.exists(outputs_dir):
        os.makedirs(outputs_dir)

    for i, input in enumerate(inputs):
        if not "mutatex" in input: continue  # skip non mutatex output files
        if input.endswith("tar.gz"): continue

        full_input_path = f"{inputs_dir}/{input}"
        full_output_path = f"{outputs_dir}/{input}.csv" # TODO how to make this .csv into something generic?

        print(f"{i}: Processing input {full_input_path}")
        preprocess_mutatex(full_input_path, full_output_path)


if __name__ == "__main__":
    from sys import argv
    if len(argv) < 4:
        print("Usage: script [single or all] [raw_results_path] [output_path]")
        print("Example 1: script single my/mutatex/results/path  my_output_file.csv")
        print("Example 2: script all my/mutatex/results  my_output_dir")
        exit(1)
    print(argv)
    if argv[1] == "single":
        print("Processing single")
        preprocess_mutatex(argv[2], argv[3])
    elif argv[1] == "all":
        print("Processing all")
        run_multi(argv[2], argv[3])
    else:
        print("Invalid arg. Must be single or all.")
        exit(1)
