
# Haemosporidian Mitochondrial Genome PacBio Pipeline (HmtG-PacBio Pipeline)

## Description
The HmtG-PacBio Pipeline is a comprehensive computational tool designed for processing Haemosporidian mitochondrial genome data from PacBio sequencing. It encompasses stages like primer detection, alignment, DNA sequence conversion to binary format, VAE model (machine learning method), clustering using DBScan, sequence correction, and local BLAST analysis.

<!-- <img src="HmtG-PacBio.jpg" alt="Flowchart Description" width="600" height="800"/> -->

## Requirements
- Python 3.x
- MAFFT v7.520
- BLAST 2.6.0
- Python Libraries: TensorFlow 2.14.0, Keras 2.14.0, BioPython, numpy, pandas, etc.

## Installation

### 1. Clone the Repository
```bash
git clone [repository-url]
```

### 2. Install Conda
Download Conda from [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/distribution) if not installed. Follow the official installation guide.

### 3. Create and Activate Conda Environment
```bash
conda create -n HmtG-PacBio python=3.9
conda activate HmtG-PacBio
```

### 4. Install Conda Packages
```bash
conda install -c bioconda mafft
conda install -c bioconda blast
```

### 5. Install Python Modules
```bash
pip3 install -r requirements.txt
```

## Usage

### Arguments
- `-rR`, `--rawReads`: Raw PacBio sequencing reads (required)
- `-pF`, `--primerF`: Forward primer (default: GATTCTCTCCACACTTCAATTCGTACTTC)<sup>1</sup>
- `-pR`, `--primerR`: Reverse primer (default: GAAGTACGAATTGAAGTGTGGAGAGAATC)<sup>1</sup>
- `-eps`, `--epsDBScan`: Epsilon for DBScan (default: 1.0)
- `-rF`, `--RemoveFiles`: Remove temporary files (default: yes)
- `-rB`, `--blastn`: Run BLASTn locally (default: yes)

### Running the Pipeline
```bash
python HmtG-PacBio-Pipeline.py -rR [path_to_raw_reads] [options]
```

## Output
Generates temporary files (aligned sequences, binary format files, cluster data) and four final output files (cluster sequences in FASTA format, training model + clustering in PNG format, genetic distances, and BLAST output files).

---

## How to cite
Pacheco, M.A., Cepeda, A.S., Miller, E.A. et al. A new long-read mitochondrial-genome protocol (PacBio HiFi) for haemosporidian parasites: a tool for population and biodiversity studies. Malar J 23, 134 (2024). https://doi.org/10.1186/s12936-024-04961-8


---
