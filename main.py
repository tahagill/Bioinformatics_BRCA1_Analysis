# This is the central controller for the entire workflow (driver of this prooject)
import os
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
from config import DATA_DIR

# Import all necessary functions
from fetch_data import fetch_sequence, fetch_protein_from_ncbi, print_sequence_details
from pairwise_align import load_sequence, shuffle_sequence, calculate_z_score, calculate_identity
from msa import trim_sequences, clean_sequence, perform_clustal_alignment
from phylogen import main as build_phylogenetic_tree
from stats import run_blast_to_xml, generate_blast_report
from validation import fetch_pfam_domains, fetch_swissprot_annotations
from visual_stats import plot_results

def run_full_analysis():
    """Executes your entire workflow with all debug prints preserved."""
    
    # 1. Fetch Data
    print("\n" + "="*50)
    print("STEP 1: Fetching Sequences from NCBI")
    print("="*50 + "\n")
    
    # DNA Sequences
    species_data = {
        "human": "NM_007294.3",
        "mouse": "NM_009764.2",
        "chimpanzee": "XM_016933586.1"
    }
    for species, accession in species_data.items():
        filename = f"{species}_brca1.fasta"
        fetch_sequence(species, accession, filename)
    
    # Protein Sequences
    protein_data = {
        "human": "NP_009225.1",
        "mouse": "NP_033894.1",
        "chimpanzee": "XP_016789075.1"
    }
    for species, protein_id in protein_data.items():
        filename = f"{species}_brca1_protein.fasta"
        fetch_protein_from_ncbi(protein_id, filename)
    
    
    print("\n=== Sequence Details ===")
    for species in species_data.keys():
        print(f"--- {species.capitalize()} DNA ---")
        print_sequence_details(f"{species}_brca1.fasta")
    for species in protein_data.keys():
        print(f"--- {species.capitalize()} Protein ---")
        print_sequence_details(f"{species}_brca1_protein.fasta")

    # 2. Pairwise Allignment
    print("\n" + "="*50)
    print("STEP 2: Pairwise Alignment")
    print("="*50 + "\n")
    
    human_protein = load_sequence("human_brca1_protein.fasta")
    mouse_protein = load_sequence("mouse_brca1_protein.fasta")
    
    
    print(f"Human BRCA1 Protein Length: {len(human_protein)}")
    print(f"Human BRCA1 Protein Start: {human_protein[:50]}")
    print(f"Mouse BRCA1 Protein Length: {len(mouse_protein)}")
    print(f"Mouse BRCA1 Protein Start: {mouse_protein[:50]}")
    
    # Alignment and Z-score
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10.0
    
    alignments = aligner.align(human_protein, mouse_protein)
    best_alignment = alignments[0]
    print("\nBest Alignment:")
    print(best_alignment)
    
    identity = calculate_identity(str(best_alignment[0]), str(best_alignment[1]))
    print(f"Biologically Accurate Percent Identity: {identity:.2f}%")
    
    z_score = calculate_z_score(human_protein, mouse_protein)
    print(f"Z-score: {z_score:.2f}")

    # 3. MSA
    print("\n" + "="*50)
    print("STEP 3: Multiple Sequence Alignment")
    print("="*50 + "\n")
    
    species = ["human", "mouse", "chimpanzee"]
    proteins = []
    for sp in species:
        record = SeqIO.read(os.path.join(DATA_DIR, f"{sp}_brca1_protein.fasta"), "fasta")
        clean_seq = clean_sequence(record.seq)
        proteins.append(SeqRecord(Seq(clean_seq), id=record.id, description=record.description))
    
    trimmed = trim_sequences(proteins, 1812)
    print("Trimmed Sequences:")
    for seq in trimmed:
        print(f"{seq.id}: Length {len(seq.seq)}, First 60 AA: {seq.seq[:60]}...")
    
    input_fasta = os.path.join(DATA_DIR, "input_sequences.fasta")
    with open(input_fasta, "w") as f:
        for seq in trimmed:
            f.write(f">{seq.id}\n{seq.seq}\n")
    
    output_fasta = os.path.join(DATA_DIR, "aligned_sequences.fasta")
    perform_clustal_alignment(input_fasta, output_fasta)

    # 4. Phylogenetics
    print("\n" + "="*50)
    print("STEP 4: Phylogenetic Tree")
    print("="*50 + "\n")
    build_phylogenetic_tree()  

    # 5. BLAST Analysis
    print("\n" + "="*50)
    print("STEP 5: BLAST Analysis")
    print("="*50 + "\n")
    
    human_protein_seq = SeqIO.read(os.path.join(DATA_DIR, "human_brca1_protein.fasta"), "fasta").seq
    blast_results = os.path.join(DATA_DIR, "blast_results.xml")
    report_path = os.path.join(DATA_DIR, "blast_report.txt")
    
    if not os.path.exists(blast_results):
        run_blast_to_xml(str(human_protein_seq), blast_results)
    
    generate_blast_report(blast_results, report_path)
    print(f"BLAST Report Generated: {report_path}")

    # 6. Validation
    print("\n" + "="*50)
    print("STEP 6: Domain Validation")
    print("="*50 + "\n")
    fetch_pfam_domains("P38398")
    fetch_swissprot_annotations("P38398")

    # 7. Visualisation
    print("\n" + "="*50)
    print("STEP 7: Visualization")
    print("="*50 + "\n")
    plot_results()

if __name__ == "__main__":
    # Verify NCBI email is set
    if Entrez.email == "your.student@university.edu":
        print("ERROR: Update your email in config.py!")
    else:
        run_full_analysis()
        print("\n=== ANALYSIS COMPLETE ===")
        print(f"All results saved to: {os.path.abspath(DATA_DIR)}")
