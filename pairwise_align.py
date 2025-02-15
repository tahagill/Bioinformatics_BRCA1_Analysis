from Bio import Align
from Bio.Align import substitution_matrices
from Bio import SeqIO
import os
import numpy as np
import random  
from config import DATA_DIR

# Load sequences from FASTA files
def load_sequence(filename):
    filepath = os.path.join(DATA_DIR, filename)
    record = SeqIO.read(filepath, "fasta")
    return str(record.seq)

# Initialize the aligner
aligner = Align.PairwiseAligner()
aligner.mode = "global"  # Global alignment
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -10.0  # More stringent gap penalty
aligner.extend_gap_score = -0.5  # Gap extension penalty

# CUSTOM SHUFFLE FUNCTION (PRESERVES AMINO ACID COMPOSITION)
def shuffle_sequence(seq):
    seq_list = list(seq)
    random.shuffle(seq_list)
    return ''.join(seq_list)

# Function to calculate Z-score
def calculate_z_score(query_seq, subject_seq, num_shuffles=1000):

    np.random.seed(42)  # Set random seed for reproducibility
    random.seed(42)  # Seed Python's random for custom shuffle

    # Calculate the original alignment score
    alignments = aligner.align(query_seq, subject_seq)
    original_score = alignments.score

    # Generate shuffled sequences and calculate their scores
    shuffled_scores = []
    for _ in range(num_shuffles):
        shuffled_seq = shuffle_sequence(subject_seq)  # Use custom shuffle
        shuffled_alignments = aligner.align(query_seq, shuffled_seq)
        shuffled_score = shuffled_alignments.score
        shuffled_scores.append(shuffled_score)
    
    # Calculate mean and standard deviation of shuffled scores
    mean = np.mean(shuffled_scores)
    std = np.std(shuffled_scores)

    # Print debugging information
    print(f"Original Score: {original_score}")
    print(f"Mean of Shuffled Scores: {mean}")
    print(f"Standard Deviation of Shuffled Scores: {std}")

    # Handle division by zero
    if std == 0:
        return float('inf')  # Indicates perfect alignment
    
    # Calculate Z-score
    z_score = (original_score - mean) / std
    return z_score

# Load and align sequences
human_protein = load_sequence("human_brca1_protein.fasta")
mouse_protein = load_sequence("mouse_brca1_protein.fasta")

# Print sequence details for debugging
print(f"Human BRCA1 Protein Length: {len(human_protein)}")
print(f"Human BRCA1 Protein Start: {human_protein[:50]}")
print(f"Mouse BRCA1 Protein Length: {len(mouse_protein)}")
print(f"Mouse BRCA1 Protein Start: {mouse_protein[:50]}")

# Perform alignment
alignments = aligner.align(human_protein, mouse_protein)

# Get the first alignment
best_alignment = alignments[0]

# Print the best alignment for debugging
print("Best Alignment:")
print(best_alignment)

# Extract aligned sequences (with gaps) directly as strings
aligned_human = str(best_alignment[0])
aligned_mouse = str(best_alignment[1])

def calculate_identity(aligned_seq1, aligned_seq2):
    matches = 0
    total = 0
    for char1, char2 in zip(aligned_seq1, aligned_seq2):
        # Only compare positions where neither sequence has a gap
        if char1 != '-' and char2 != '-':
            total += 1
            if char1 == char2:
                matches += 1
    
    if total == 0:
        return 0
    return (matches / total) * 100

identity = calculate_identity(aligned_human, aligned_mouse)
print(f"Biologically Accurate Percent Identity: {identity:.2f}%")

# Calculate Z-score
z_score = calculate_z_score(human_protein, mouse_protein)
print(f"Z-score for human-mouse alignment: {z_score:.2f}")