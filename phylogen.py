import os
import matplotlib.pyplot as plt
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator
from config import DATA_DIR

def main():
    # Load aligned sequences
    aligned_file = os.path.join(DATA_DIR, "aligned_sequences.fasta")
    alignment = AlignIO.read(aligned_file, "fasta")

    # Calculate distance matrix
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(alignment)

    # Build tree
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    for clade in tree.find_clades():
        if clade.name:
            clade.name = clade.name.replace("NP_009225.1", "Human") \
                                   .replace("NP_033894.1", "Mouse") \
                                   .replace("XP_016789075.1", "Chimpanzee")

    # Ensure non-zero branch lengths (critical for visualization)
    for clade in tree.find_clades():
        if clade.branch_length is None or clade.branch_length < 0.0001:
            clade.branch_length = 0.0001  # Set minimal branch length

    try:
        plt.figure(figsize=(12, 8))
        ax = plt.subplot(1, 1, 1)
        
        Phylo.draw(
            tree,
            axes=ax,
            branch_labels=lambda c: f"{c.branch_length:.4f}" if c.branch_length else "",
            label_colors=lambda _: "black",  # Force black labels
            do_show=False
        )
        
        # Adjust layout
        plt.title("BRCA1 Phylogenetic Tree", fontsize=14)
        plt.tight_layout() 
        plt.show()
        
    except Exception as e:
        print(f"Visualization failed: {e}. Using ASCII fallback:")
        Phylo.draw_ascii(tree)
