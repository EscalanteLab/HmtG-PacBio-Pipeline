from Bio import SeqIO
from re import search
import random


def primerDetection(input_file, primerF, primerR):
    # Read all sequences from a fasta file and store them in a dictionary with their IDs as keys and sequences as values
    seq_dict = {rec.id: rec.seq for rec in SeqIO.parse(input_file, "fasta")}
    outf = str(input_file).split(".")[0]
    # If the dictionary contains more than 5000 sequences, randomly select 5000 sequences
    if len(seq_dict) > 5000:
        # seq_dict = dict(random.sample(seq_dict.items(), 5000))
        seq_dict = dict(random.sample(list(seq_dict.items()), 5000))

    # Open the output files using "with open", this ensures files are closed properly even if an error occurs
    with open(outf + "_nolong.fasta", "a") as f, open(outf + "_long.fasta", "a") as g:
        for id, seq in seq_dict.items():  # Iterate over each sequence in the dictionary
            seq_str = str(seq)  # Convert sequence to string
            rev_com = (
                seq.reverse_complement()
            )  # Calculate reverse complement of sequence
            # Determine the output file based on sequence length
            outfile = f if len(seq_str) <= 6000 else g

            # If the forward primer is found in the sequence
            if search(primerF, seq_str):
                # Write the sequence ID and sequence to the appropriate output file
                print(">" + id, file=outfile)
                print(seq_str, file=outfile)
            # If the reverse primer is found in the sequence
            elif search(primerR, seq_str):
                # Write the sequence ID and reverse complement sequence to the appropriate output file
                print(">" + id, file=outfile)
                print(str(rev_com), file=outfile)
