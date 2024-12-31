# utilities
import subprocess
from Bio import SeqIO
import os
import glob
import numpy as np
import pandas as pd
import shutil

# Define the database prefix (adjust this based on your actual database path)
db_prefix = 'sprot/sprot'

def run_PSIBLAST(input_file, out_dir):
    """
    This function runs PSI-BLAST for each sequence in the input FASTA file and saves the PSSM profiles.
    :param input_file: Path to the input FASTA file.
    :param out_dir: Directory where PSSM profiles will be saved.
    """
    os.makedirs(out_dir, exist_ok=True)

    # Open the log file to store command output and errors
    with open('log_file.txt', 'w') as log:
        try:
            # Parse each sequence in the FASTA file
            for record in SeqIO.parse(input_file, 'fasta'):
                name = record.id  # Protein name or identifier
                sequence = str(record.seq)  # Protein sequence

                # Write sequence to a temporary file in FASTA format
                temp_file = 'temp.txt'
                with open(temp_file, 'w') as temp:
                    temp.write(f'>{name}\n{sequence}\n')

                # Define the PSIBLAST command
                command = [
                    'psiblast',
                    '-query', temp_file,
                    '-db', db_prefix,
                    '-num_iterations', '3',
                    '-out_ascii_pssm', os.path.join(out_dir, f"{name}.pssm")
                ]

                # Run the command and capture output in the log file
                try:
                    print(f"Running PSI-BLAST for: {name}")
                    subprocess.run(command, stdout=log, stderr=log, text=True, check=True)
                except subprocess.CalledProcessError as e:
                    print(f"Error running PSI-BLAST for {name}. See log file for details.")
                    log.write(f"\nError occurred while processing {name}:\n{e}\n")

                # Remove the temporary file after use
                os.remove(temp_file)

        except Exception as e:
            # Handle unexpected exceptions
            print(f"An unexpected error occurred: {e}")
            log.write(f"\nUnexpected error: {e}\n")


def load_pssm_file(filename):
    data = []
    with open(filename, 'r') as file:
        for line in file:
            # Split line by whitespace
            parts = line.strip().split()
            # Check if the line has exactly 44 columns (adjust as per actual data structure)
            if len(parts) == 44:
                # Extract the relevant columns (e.g., columns 2 to 21 for a 20-column PSSM matrix)
                row_data = list(map(float, parts[2:22]))  # Convert to float
                data.append(row_data)
    # Convert the list to a numpy array
    return np.array(data)


def bigram_pssm(pssm_dir):

    #search for all .pssm files in the output directory
    pssm_files = glob.glob(pssm_dir+"*.pssm")
    # Step 1: Load PSSM files into a list of numpy arrays
    pssm_data = []
    for file in pssm_files:
        # Load the data
        data = load_pssm_file(file)
        pssm_data.append(data)

    # Step 2: Extract and process matrices for each protein sequence
    processed_data = []
    for data in pssm_data:
        # Exclude the last 5 rows and keep only the first 20 columns
        matrix = data[:-5, :20]  # M-5 rows and first 20 columns

        # Note: No normalization is applied as per the MATLAB code
        processed_data.append(matrix)
    # Step 3: Calculate Bi-PSSM features
    bi_pssm_features = []
    
    for matrix in processed_data:
        M, N = matrix.shape
        xx = np.zeros((M - 1, 20 * 20))  # Prepare matrix for storing pairwise products
        p = 0  # Initialize feature index

        # Nested loops to calculate pairwise products for consecutive rows
        for m in range(20):
            for n in range(20):
                # Compute product of (i, m) and (i+1, n) elements across rows
                xx[:, p] = matrix[:-1, m] * matrix[1:, n]
                p += 1

        # Sum over rows for each pairwise feature and store in `B`
        bi_pssm_features.append(np.sum(xx, axis=0))

    # headers
    amino_acids = 'ARNDCQEGHILKMFPSTWYV'
    header = [f'BiG_PSSM_{aa1}{aa2}' for aa1 in amino_acids for aa2 in amino_acids]

    # Convert results to a numpy array and then to a DataFrame
    bi_pssm_features = np.array(bi_pssm_features)
    bi_pssm_df = pd.DataFrame(bi_pssm_features, columns=header)

    bi_pssm_df.insert(0, 'Protein_ID', [file[15:-5] for file in pssm_files])

    # Step 4: Save the results to a file
    bi_pssm_df.to_excel("bipssmresults.xlsx", index=False, header=True)




def main():
    # Input FASTA file and output directory for PSSM profiles
    input_file = 'test_fasta.txt'
    out_dir = 'output_BiGpssm/'

    # Check if input file exists
    if not os.path.isfile(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        return

    # Run PSI-BLAST and process PSSM profiles
    run_PSIBLAST(input_file, out_dir)
    bigram_pssm(out_dir)

    # Clean up the output directory
    try:
        shutil.rmtree(out_dir, ignore_errors=True)
    except Exception as e:
        print(f"Error cleaning up directory '{out_dir}': {e}")


if __name__ == "__main__":
    main()
