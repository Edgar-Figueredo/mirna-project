#  Integrative miRNA Annotation Pipeline  
### Cross-database harmonization of miRNA gene annotations and sequences
---
##  Overview
This project implements a reproducible bioinformatics pipeline to integrate and compare microRNA (miRNA) annotations from multiple curated databases:
- GENCODE → gene-centric genomic annotation  
- miRBase → sequence-centric mature miRNAs  
- miRGeneDB → high-confidence curated miRNA genes  
The pipeline harmonizes naming conventions, extracts genomic sequences, and performs cross-database comparisons to identify consensus and database-specific miRNAs.
---
##  Objectives
- Standardize miRNA nomenclature across databases  
- Extract genomic coordinates and sequences of miRNAs  
- Compare gene-level annotations between datasets  
- Identify high-confidence miRNAs supported across resources  
- Generate quality control (QC) reports and summary statistics  
---
##  Biological Context
miRNA annotation varies significantly across databases:
- GENCODE provides genomic loci but lacks mature sequences  
- miRBase provides mature sequences but includes varying confidence levels  
- miRGeneDB focuses on evolutionarily conserved, high-confidence miRNAs  
This pipeline integrates these complementary layers to produce a harmonized and biologically meaningful miRNA dataset.
---
##  Pipeline Structure
The workflow consists of the following steps:
1. Data acquisition
   - Download GENCODE, miRBase, and miRGeneDB data
2. Parsing
   - Extract miRNA annotations from GENCODE (GTF)
   - Extract human miRNAs from miRBase (FASTA)
   - Parse curated miRNAs from miRGeneDB
3. Coordinate processing
   - Convert GENCODE annotations to BED format
   - Remove duplicated genomic loci
4. Sequence extraction
   - Retrieve genomic sequences using BEDTOOLS
5. Normalization
   - Harmonize miRNA naming (e.g., hsa-miR-21-5p → MIR21)
6. Cross-database comparison
   - Identify shared and unique miRNAs
7. Quality Control
   - Validate FASTA headers
   - Generate QC report
8. Visualization
   - Venn diagram of dataset overlap
---
##  Project Structure
<img width="352" height="243" alt="image" src="https://github.com/user-attachments/assets/0932652f-8b38-425c-abe9-102325affbf4" />

---
## 🚀 Installation

### 1. Clone repository and enter the project
```bash
git clone https://github.com/Edgar-Figueredo/mirna-project.git
cd mirna-project
```

---

### 2. Setup environment

> **The setup scripts check for Python automatically. If Python is not installed, it will be downloaded and installed for you.**

---

#### 🍎 Mac / Linux
```bash
bash setup.sh
```

---

#### 🪟 Windows (CMD) — recommended for most users

Simply double-click `setup.bat`, or run it from the Command Prompt:

```cmd
setup.bat
```

> If Python is not installed, the script will **download and install Python 3.11 automatically**.  
> After installation completes, **close the window and run `setup.bat` again** so the PATH is refreshed.

---

#### 🪟 Windows (PowerShell)

If you prefer PowerShell, you may need to allow script execution first (one-time):

```powershell
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
```

Then run:

```powershell
.\setup.ps1
```

> Same as above: if Python is missing, it will be installed automatically.  
> Restart PowerShell and run `.\setup.ps1` again if prompted.

---

### 3. Activate environment (for future sessions)

Once setup is done, the environment is activated automatically. To reactivate it later:

#### 🍎 Mac / Linux
```bash
source miRNA_env/bin/activate
```

#### 🪟 Windows (PowerShell)
```powershell
.\miRNA_env\Scripts\Activate.ps1
```

#### 🪟 Windows (CMD)
```cmd
miRNA_env\Scripts\activate.bat
```

---

### 4. Run pipeline
```bash
python scripts/mirna_super_pipeline.py
```

---

##  Usage

Run the full pipeline:

```bash
python scripts/mirna_super_pipeline.py
```

---

##  Outputs

All results are stored in `data/processed/`

Key outputs include:

| File | Description |
|------|------------|
| `gencode_mirna.txt` | Raw miRNA gene names |
| `mirna_coordinates_unique.bed` | Genomic miRNA loci |
| `mirna_sequences.fa` | Extracted genomic sequences |
| `gencode_norm.txt` | Normalized GENCODE miRNAs |
| `mirbase_norm.txt` | Normalized miRBase miRNAs |
| `final_overlap.txt` | Shared and unique miRNAs |
| `qc_report.txt` | Quality control summary |
| `venn_mirna.png` | Overlap visualization |

---

##  Quality Control

The pipeline automatically evaluates:
- Number of extracted sequences  
- FASTA header normalization  
- Cross-database agreement  

Example:
```
FASTA sequences: 1854   Bad headers: 0   Database overlap: HIGH
```

---

##  Key Insights
- miRNA annotation is database-dependent and inconsistent  
- Normalization is essential for meaningful comparison  
- miRGeneDB provides a high-confidence subset  
- Overlap between GENCODE and miRBase is substantial but incomplete  

---

##  Limitations
- No miRNA-target interaction data included  
- No expression or functional annotation  
- miRBase includes entries with variable confidence  

---

##  Future Work
- Integrate target databases:
  - TargetScan  
  - miRTarBase  
  - TarBase  
- Incorporate expression datasets (RNA-seq / small RNA-seq)  
- Add confidence scoring (multi-database support)  
- Extend to other species  

---

##  Author

Edgar Alfonso Figueredo Carmona  
Microbiologist & Bioanalysis, MSc student: bioinformatics  
Universidad de Antioquia  

---

##  License

This project is intended for academic and research purposes.
