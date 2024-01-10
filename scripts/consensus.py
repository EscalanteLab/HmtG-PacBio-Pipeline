from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import os


def cons(list_files):
    for input_file in list_files:
        nameOut = str(input_file).split("_mafft")[0]
        outf = str(input_file).split("_mafft")[0] + "_consensus.fasta"
        alignment = AlignIO.read(input_file, "fasta")  # Read the alignment file
        consensus_seq = ""  # Initialize consensus sequence

        for i in range(len(alignment[0])):  # Loop over the alignment length
            nucleotide_counts = {}  # Dictionary to count occurrence of each nucleotide
            gap_counts = 0  # Counter for gaps

            for record in alignment:  # For each sequence in the alignment
                nucleotide = record.seq[i]  # Get the nucleotide at the current position
                nucleotide_counts[nucleotide] = (
                    nucleotide_counts.get(nucleotide, 0) + 1
                )  # Count the nucleotide

                # Count the gaps
                if nucleotide == "-":
                    gap_counts += 1

            # Continue to the next position if more than 50% are gaps
            if gap_counts / len(alignment) > 0.5:
                continue

            # Sort the nucleotides by counts and get the most common nucleotide
            sorted_nucleotides = sorted(
                nucleotide_counts.items(), key=lambda x: x[1], reverse=True
            )

            # Check if the most common nucleotide appears in more than 95% of the sequences
            if sorted_nucleotides[0][1] / len(alignment) > 0.85:
                consensus_seq += sorted_nucleotides[0][
                    0
                ]  # Add the most common nucleotide to the consensus sequence
            else:
                consensus_seq += "N"  # Add 'X' if there is no consensus

        # Write the consensus sequence to a new fasta file
        consensus_seq = Seq(consensus_seq)
        consensus_record = SeqRecord(consensus_seq, id=nameOut, description="")
        SeqIO.write(consensus_record, outf, "fasta")
