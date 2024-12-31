# utilities
import subprocess
from Bio import SeqIO
import os

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

def main():
    # Input FASTA file and output directory for PSSM profiles
    input_file = ''         #enter the input file directory
    out_dir = ''   #enter the output directory

    run_PSIBLAST(input_file, out_dir)

if __name__ == "__main__":
    main()
