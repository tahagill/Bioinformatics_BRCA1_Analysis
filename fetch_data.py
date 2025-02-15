import os
import time
from Bio import Entrez, SeqIO
from config import DATA_DIR, Entrez


os.makedirs(DATA_DIR, exist_ok=True)

# List of species and their BRCA1 accession numbers (DNA)
species_data = {
    "human": "NM_007294.3",      # Human BRCA1 mRNA
    "mouse": "NM_009764.2",      # Mouse BRCA1 mRNA
    "chimpanzee": "XM_016933586.1"  # Chimpanzee BRCA1 mRNA
}

# Function to fetch and save sequences
def fetch_sequence(species, accession, filename):
    try:
        # Fetch the sequence from GenBank
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        
        # Save the sequence in FASTA format in the DATA_DIR
        filepath = os.path.join(DATA_DIR, filename)
        SeqIO.write(record, filepath, "fasta")
        print(f"Saved {species} BRCA1 sequence to {filepath}")
    except Exception as e:
        print(f"Error fetching {species} sequence: {e}")

# Fetch sequences for all species
for species, accession in species_data.items():
    filename = f"{species}_brca1.fasta"
    fetch_sequence(species, accession, filename)

# List of species and their NCBI Protein IDs
protein_data = {
    "human": "NP_009225.1",      # Human BRCA1 protein
    "mouse": "NP_033894.1",      # Mouse BRCA1 protein
    "chimpanzee": "XP_016789075.1"  # Chimpanzee BRCA1 protein
}

# Function to fetch protein sequence from NCBI
def fetch_protein_from_ncbi(protein_id, filename):
    try:
        # Fetch the protein sequence from NCBI
        handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        
        # Save the sequence in FASTA format in the DATA_DIR
        filepath = os.path.join(DATA_DIR, filename)
        SeqIO.write(record, filepath, "fasta")
        print(f"Saved protein sequence {protein_id} to {filepath}")
    except Exception as e:
        print(f"Error fetching protein {protein_id}: {e}")

# Fetch protein sequences for all species
for species, protein_id in protein_data.items():
    filename = f"{species}_brca1_protein.fasta"
    fetch_protein_from_ncbi(protein_id, filename)
    time.sleep(1)  # Added a small delay between requests to avoid rate limits

# Function to print sequence details
def print_sequence_details(filename):
    filepath = os.path.join(DATA_DIR, filename)
    if os.path.exists(filepath):
        record = SeqIO.read(filepath, "fasta")
        print(f"Sequence ID: {record.id}")
        print(f"Description: {record.description}")
        print(f"Sequence Length: {len(record.seq)}")
        print(f"First 50 bases: {record.seq[:50]}\n")
    else:
        print(f"File not found: {filepath}")

# Print details for all sequences
for species in species_data.keys():
    print(f"--- {species.capitalize()} BRCA1 DNA Sequence ---")
    print_sequence_details(f"{species}_brca1.fasta")

for species in protein_data.keys():
    print(f"--- {species.capitalize()} BRCA1 Protein Sequence ---")
    print_sequence_details(f"{species}_brca1_protein.fasta")