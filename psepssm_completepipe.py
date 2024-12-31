# utilities
import subprocess
from Bio import SeqIO
import os
import glob
import numpy as np
import pandas as pd
import shutil
import csv

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


def psepssm_calc(pssm_dir, num_lags):
    """
    Calculate PsePSSM for a list of PSSM matrices and save to a CSV file.

    Parameters:
        pssm_matrices (list of numpy arrays): List of PSSM matrices (M x 20).
        protein_names (list of str): List of protein names.
        num_lags (int): Number of lag parameters (λ).
        output_file (str): Path to save the CSV file.
    """
    aavec = 'ARNDCQEGHILKMFPSTWYV'

    #if len(pssm_matrices) != len(protein_names):
        #raise ValueError("The number of proteins and PSSM matrices must match.")


    # search for all .pssm files in the output directory
    pssm_files = glob.glob(pssm_dir + "*.pssm")
    # Step 1: Load PSSM files into a list of numpy arrays
    pssm_matrices = []
    for file in pssm_files:
        # Load the data
        data = load_pssm_file(file)
        pssm_matrices.append(data)

    num_proteins = len(pssm_matrices)
    normalized_matrices = []

    #Step 2 : Normalize PSSM matrices
    for matrix in pssm_matrices:
        normalized_matrix = 1 / (1 + np.exp(-matrix))  # Logistic normalization
        normalized_matrices.append(normalized_matrix)

    # Step 3 : Calculate PSSM-AAC
    pssm_aac = np.zeros((num_proteins, 20))
    for i, matrix in enumerate(normalized_matrices):
        pssm_aac[i, :] = np.mean(matrix, axis=0)

    # Step 4 Calculate PsePSSM components
    pse_pssm_features = []
    for lag in range(1, num_lags + 1):
        lag_features = np.zeros((num_proteins, 20))
        for t, matrix in enumerate(normalized_matrices):
            num_rows, _ = matrix.shape
            for j in range(20):
                squared_differences = [
                    (matrix[i, j] - matrix[i + lag, j]) ** 2 for i in range(num_rows - lag)
                ]
                lag_features[t, j] = sum(squared_differences) / (num_rows - lag)
        pse_pssm_features.append(lag_features)

    pse_pssm_features = np.hstack(pse_pssm_features)  # Combine all lag features
    psepssm = np.hstack([pssm_aac, pse_pssm_features])  # Concatenate AAC and PsePSSM features

    # Create headers for the CSV
    amino_acids = [f"Normalized_PSSM_AAC_{aa}" for aa in aavec]  # AAC headers (1-20)
    lag_headers = [
        f"PsePSSM_l{lag}_{aa}" for lag in range(1, num_lags + 1) for aa in aavec
    ]
    headers = ["Protein_ID"] + amino_acids + lag_headers
    protein_names = [file[15:-5] for file in pssm_files]
    # Write to CSV
    with open('psepssm_results.csv', mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(headers)  # Write header
        for i, row in enumerate(psepssm):
            data = [protein_names[i]]
            data.extend(row.tolist())
            # data.extend(row.tolist())
            # print(data)
            writer.writerow(data)  # Write protein name and features

    #print(f"PsePSSM features saved to {output_file}")


def main():
    # Input FASTA file and output directory for PSSM profiles
    input_file = 'test_fasta.txt'
    out_dir = 'output_psepssm/'
    num_lags = 3

    # Check if input file exists
    if not os.path.isfile(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        return

    # Run PSI-BLAST and process PSSM profiles
    run_PSIBLAST(input_file, out_dir)
    psepssm_calc(out_dir,num_lags)

    # Clean up the output directory
    try:
        shutil.rmtree(out_dir, ignore_errors=True)
    except Exception as e:
        print(f"Error cleaning up directory '{out_dir}': {e}")

if __name__ == "__main__":
    main()
