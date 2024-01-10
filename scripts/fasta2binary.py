from Bio import SeqIO  # Importing SeqIO from BioPython for parsing sequences.
import subprocess

def nucleotide_to_binary(nucleotide):
    nucleotide = nucleotide.upper()
    nucleotide_dict = {
        "A": "1,0,0,0",
        "C": "0,1,0,0",
        "G": "0,0,1,0",
        "T": "0,0,0,1",
        "U": "0,1,1,0",  # Uracil
        "R": "1,1,0,0",  # A o G
        "Y": "0,1,1,1",  # C o T
        "S": "0,1,0,1",  # G o C
        "W": "1,0,1,0",  # A o T
        "M": "1,0,1,1",  # A o C
        "B": "0,0,1,1",  # C o G o T
        "D": "1,1,0,1",  # A o G o T
        "H": "1,0,0,1",  # A o C o T
        "V": "1,0,1,1",  # A o C o G
        "N": "0,0,0,0",  # Any base
        "-": "0,0,0,0"   # Gap
    }

    # Return the binary representation if the nucleotide is found in the dictionary
    if nucleotide in nucleotide_dict:
        return nucleotide_dict[nucleotide]
    else:
        print("Invalid nucleotide provided")

# This function converts a multifasta file to a binary sequence file.
def fasta2bin(input_file, dir_out):
    outf0 = str(input_file).split("/")[-1]
    nameFile = str(outf0).split("_")[0]
    output_file = str(outf0).split("_")[0] + "_bin.txt"


    with open(input_file, 'r') as fasta_file, open(output_file, 'w') as output:
        for record in SeqIO.parse(fasta_file, 'fasta'):  # Iterate over each record in the fasta file.
            header = record.id
            sequence = record.seq
            # Convert each nucleotide in the sequence to its binary representation and join them by spaces.
            binary_sequence = ' '.join([nucleotide_to_binary(nt) for nt in sequence])
            output.write(header + " " + nameFile + " " + binary_sequence + "\n")  # Write the header and binary sequence to the output file.
    mv_command = ["mv", output_file, dir_out]
    subprocess.run(mv_command)


