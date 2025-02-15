from Bio import SeqIO
from Bio.Seq import Seq  
from Bio.SeqRecord import SeqRecord
import os
import subprocess
from config import DATA_DIR
# Define the target length for trimming (1812 for the mouse sequence)
TARGET_LENGTH = 1812

def trim_sequences(sequences, target_length):
    trimmed_sequences = []
    for seq in sequences:
        # Trim each sequence to the target length
        trimmed_seq = seq.seq[:target_length]
        trimmed_sequences.append(SeqRecord(trimmed_seq, id=seq.id, description=seq.description))
    return trimmed_sequences

def clean_sequence(seq):
    # Convert to string and uppercase
    seq_str = str(seq).upper()
    # Keep only valid amino acid characters
    valid_chars = set("ACDEFGHIKLMNPQRSTVWY-")
    return ''.join(c for c in seq_str if c in valid_chars)

# Load and clean sequences
species = ["human", "mouse", "chimpanzee"]
proteins = []
for sp in species:
    filename = f"{sp}_brca1_protein.fasta"
    filepath = os.path.join(DATA_DIR, filename)
    record = SeqIO.read(filepath, "fasta")
    
    # Clean the sequence
    clean_seq = clean_sequence(record.seq)
    # Create new record with cleaned sequence
    new_record = SeqRecord(
        Seq(clean_seq),
        id=record.id,
        description=record.description
    )
    proteins.append(new_record)

# Trim the sequences to the target length (1812 for mouse)
trimmed_sequences = trim_sequences(proteins, TARGET_LENGTH)

# Print trimmed sequence information
print("Loaded and Trimmed sequences:")
for seq in trimmed_sequences:
    print(f"\n{seq.id}")
    print(f"Length: {len(seq.seq)}")
    print(f"First 60 aa: {seq.seq[:60]}")

# Save the trimmed sequences to new files in DATA_DIR
for seq in trimmed_sequences:
    trimmed_filepath = os.path.join(DATA_DIR, f"trimmed_{seq.id}.fasta")
    with open(trimmed_filepath, "w") as handle:
        SeqIO.write(seq, handle, "fasta")
    print(f"Trimmed sequence for {seq.id} saved to {trimmed_filepath}")

# Perform MSA using Clustal Omega
def perform_clustal_alignment(input_fasta, output_fasta):
    try:
        # Running Clustal Omega with subprocess 
        command = [
            "clustalo",  # The Clustal Omega executable
            "-i", input_fasta,  # Input FASTA file
            "-o", output_fasta,  # Output FASTA file for aligned sequences
            "--force"  # Force overwriting the output file if it exists
        ]
        
        # Call Clustal Omega via subprocess
        subprocess.run(command, check=True)
        print(f"Alignment complete. Results saved to {output_fasta}")
    except subprocess.CalledProcessError as e:
        print(f"Error running Clustal Omega: {e}")

# Save the trimmed sequences to an input file for Clustal Omega
input_fasta = os.path.join(DATA_DIR, "input_sequences.fasta")
with open(input_fasta, "w") as handle:
    for seq in trimmed_sequences:
        handle.write(f">{seq.id}\n{str(seq.seq)}\n")
print(f"Trimmed sequences saved to {input_fasta}")

# Define the output file for Clustal Omega alignment
output_fasta = os.path.join(DATA_DIR, "aligned_sequences.fasta")

# Perform the alignment using Clustal Omega
perform_clustal_alignment(input_fasta, output_fasta)

# Print the aligned sequences from the output file
try:
    with open(output_fasta, "r") as handle:
        aligned_sequences = handle.read()
    print("\nAligned Sequences:")
    print(aligned_sequences)
except Exception as e:
    print(f"Error reading the aligned sequences: {e}")
