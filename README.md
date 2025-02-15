# BRCA1 Evolutionary Insights

## 📌 Description
This project analyzes **BRCA1 protein sequences** from humans, mice, and chimpanzees using **Biopython**, **Clustal Omega**, and phylogenetic tools. It performs:
- **Pairwise alignment**
- **Multiple sequence alignment (MSA)**
- **Phylogenetic tree construction**
- **Z-score calculations**
- **BLAST searches**
- **Domain validation using UniProt**
- **BLAST Report sheet in an awesome format**

The goal is to study **evolutionary conservation** and identify functionally important regions in the BRCA1 protein.

---

## 🛠 Tools & Libraries Used
- **Biopython** → Sequence alignment, BLAST, phylogenetic tree construction
- **Clustal Omega** → Multiple sequence alignment
- **Matplotlib / Seaborn** → Data visualization
- **NCBI Entrez** → Fetching sequences
- **UniProt & Swissprot** → Domain validation & functional annotations

---

## 🎯 Project Objectives
1. **Fetch BRCA1 protein sequences** for humans, mice, and chimpanzees.
2. **Perform pairwise alignment** to calculate **percent identity** and **Z-scores**.
3. **Conduct multiple sequence alignment (MSA)** using Clustal Omega.
4. **Build a phylogenetic tree** to visualize evolutionary relationships.
5. **Validate conserved domains** using UniProt annotations.
6. **Visualize results** using Matplotlib & Seaborn.

---

## 🚀 Installation & Setup

### 1️⃣ Clone the Repository
```bash
git clone https://github.com/tahagill/Bioinformatics_BRCA1_Analysis.git
cd BRCA1_Evolutionary_Insights
```

### 2️⃣ Set Up a Virtual Environment
```bash
python -m venv venv
```
#### Activate the Virtual Environment:
- **Windows**
  ```bash
  venv\Scripts\activate
  ```
- **macOS/Linux**
  ```bash
  source venv/bin/activate
  ```

### 3️⃣ Install Dependencies
```bash
pip install -r requirements.txt
```

### 4️⃣ Run the Analysis
```bash
python main.py
```

---

## 📂 Project Structure
```
BRCA1_Evolutionary_Insights/
├── data/                   # Auto-created folder for output files
├── src/                    # Source code
│   ├── config.py           # Configuration file (email, paths)
│   ├── fetch_data.py       # Fetches sequences from NCBI
│   ├── pairwise_align.py   # Performs pairwise alignment and Z-score calculation
│   ├── msa.py              # Handles multiple sequence alignment
│   ├── phylogen.py         # Builds phylogenetic trees
│   ├── stats.py            # Runs BLAST and generates reports
│   ├── validation.py       # Validates domains using UniProt
│   └── visual_stats.py     # Visualizes results
├── .gitignore              # Ignores data/, venv/, and cache files
├── main.py                 # Entry point for the analysis
├── requirements.txt        # List of dependencies
└── README.md               # Project documentation
```

---

## 🔧 Configuration (`config.py`)
The `config.py` file contains the following configurations:
- **NCBI Email** → Required for API access.
- **DATA_DIR** → Automatically creates a data folder for output files.

### Editing `config.py`
Modify the following lines if you want to use your own NCBI email or API key:
```python
Entrez.email = "your.email@university.edu"  # Replace with your email
Entrez.api_key = "your_api_key_here"       # Optional: Add your API key
```

---

## 🔬 Key Analyses
### **1️⃣ Pairwise Alignment**
- Compares human BRCA1 with mouse and chimpanzee sequences.
- Calculates **percent identity** and **Z-scores** to assess alignment significance.

### **2️⃣ Multiple Sequence Alignment (MSA)**
- Aligns BRCA1 sequences from multiple species using Clustal Omega.
- Trims sequences for equal length analysis.

### **3️⃣ Phylogenetic Tree**
- Constructs a **Neighbor-Joining (NJ) tree** based on sequence identity.
- Visualizes evolutionary relationships between species.

### **4️⃣ BLAST Analysis**
- Runs **BLASTp** on human BRCA1 to identify homologous sequences.
- Generates detailed alignment reports and domain annotations.

### **5️⃣ Domain Validation**
- Fetches domain annotations from **UniProt**.
- Identifies key conserved regions like **RING** and **BRCT domains**.

### **6️⃣ Visualization**
- Generates phylogenetic trees, BLAST results, and alignment plots using **Matplotlib** and **Seaborn**.

---

## ✅ Requirements
- **Python 3.8+**
- **Biopython**
- **NumPy**
- **Matplotlib**
- **Seaborn**
- **Clustal Omega** _(must be installed separately and added to PATH)_

---

## ⚠️ Notes
- The `data/` folder is auto-created to store output files (e.g., FASTA files, alignment results, trees).
- **Clustal Omega** must be installed and added to your system **PATH**.
  - If you encounter issues, ensure Clustal Omega is correctly installed and accessible.

---

## 📢 Contributing
Feel free to contribute by submitting **pull requests** or **opening issues**. Suggestions for improving the analysis and code are welcome!

---

