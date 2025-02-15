import matplotlib.pyplot as plt
import seaborn as sns

def plot_results():
    """Plots BLAST results: % Identity vs. E-value."""
    blast_results = [
        {"species": "Homo sapiens", "identity": 98, "evalue": 1e-50, "name": "BRCA1"},
        {"species": "Mus musculus", "identity": 85, "evalue": 1e-40, "name": "BRCA1"},
        {"species": "Rattus norvegicus", "identity": 75, "evalue": 1e-30, "name": "BRCA1"},
        {"species": "Canis lupus", "identity": 65, "evalue": 1e-20, "name": "BRCA1"},
        {"species": "Bos taurus", "identity": 55, "evalue": 1e-10, "name": "BRCA1"}
    ]

    # Extracting Data
    species_list = list(set([entry["species"] for entry in blast_results]))
    colors = sns.color_palette("tab10", len(species_list))  # Unique colors for each species
    species_color_map = {species: colors[i] for i, species in enumerate(species_list)}

    fig, ax = plt.subplots(figsize=(12, 6))

    for entry in blast_results:
        species = entry["species"]
        color = species_color_map[species]
        ax.scatter(entry["evalue"], entry["identity"], color=color, label=species, s=100, edgecolors='black', alpha=0.7)
        
        # Shorten names for clarity
        label = f"{entry['name']} ({species.split()[0]})"
        ax.text(entry["evalue"], entry["identity"], label, fontsize=10, ha='right', rotation=15)

    ax.set_xscale("log")  # Log scale for better visualization
    ax.set_xlabel("E-value (log scale)")
    ax.set_ylabel("% Identity")
    ax.set_title("BLAST Results: % Identity vs. E-value")

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))  # Remove duplicates
    ax.legend(by_label.values(), by_label.keys(), loc="upper left", bbox_to_anchor=(1, 1))

    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.tight_layout()
    plt.show()

