
# Haemosporidian Mitochondrial Genome PacBio Pipeline (HmtG-PacBio Pipeline)

## Description
The HmtG-PacBio Pipeline is a comprehensive computational tool designed for processing Haemosporidian mitochondrial genome data from PacBio sequencing. It includes various stages such as primer detection, alignment, DNA sequence conversion to binary format, VAE model (machine learning method), clustering using DBScan, correction sequence, and local BLAST analysis.

<img src="HmtG-PacBio.jpg" alt="Flowchart Description" width="600" height="800"/>

## Requirements
- Python 3.x
- MAFFT v7.520
- BLAST 2.6.0
- Various Python Libraries: `TensorFlow 2.14.0`, `Keras 2.14.0`, `BioPython`, `numpy`, `pandas`, etc.

## Installation
1. Clone the repository:
   ```
   git clone [repository-url]
   ```
2. This project utilizes Conda, a package and environment management system, which simplifies the installation and management of software packages. Below are the steps to install all necessary requirements to run this project, including Conda packages like MAFFT and BLAST, as well as Python modules.

i) Installing Conda
If you don't have Conda installed, you can download it from Miniconda or Anaconda. Follow the instructions on the official website to install it on your operating system.

ii) Creating and Activating a Conda Environment
Once Conda is installed, create a specific environment for this project. This ensures that the dependencies for this project do not interfere with those of other projects. Use the following command to create an environment:

```
conda create -n HmtG-PacBio python=3.9
```

To activate the environment, use:

```
conda activate HmtG-PacBio
```

iii) Installing Conda Packages
Within the activated Conda environment, install necessary Conda packages like MAFFT and BLAST using the following commands:

```
conda install -c bioconda mafft
conda install -c bioconda blast
```

iv) Installing Python Modules

Next, install the required Python modules using the requirements.txt file provided in this repository:

```
pip3 install -r requirements.txt
```

## Usage

### Arguments
- `-rR`, `--rawReads`: Raw reads from PacBio sequencing (required)
- `-pF`, `--primerF`: Forward primer 5'-3' (default: GATTCTCTCCACACTTCAATTCGTACTTC)<sup>1</sup>
- `-pR`, `--primerR`: Reverse primer 3'-5' (default: GAAGTACGAATTGAAGTGTGGAGAGAATC)<sup>1</sup>
- `-eps`, `--epsDBScan`: Epsilon value for DBScan clustering (default: 1.0)
- `-rF`, `--RemoveFiles`: Option to remove temporary files (default: yes)
- `-rB`, `--blastn`: Option to run BLASTn locally (default: yes)

```
References:

1. Haemosporidian Mitochondrial Genome PacBio Pipeline (HmtG-PacBio Pipeline) Reference
```


### Running the Pipeline
Run the pipeline with the required arguments:
```
python HmtG-PacBio-Pipeline.py -rR [path_to_raw_reads] [options]
```
You can also specify other optional arguments as needed.

## Output
The pipeline generates temporay files at different stages (including aligned sequences, binary format files, cluster data, etc) and four final output files (cluster sequences in fasta format, training model + clustering in a png format, genetics distances and blast output files).


