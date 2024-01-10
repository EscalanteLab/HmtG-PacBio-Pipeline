from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from collections import Counter


def countHaplotypes(alignment_file, output_file):
    # Read the alignment from the provided file in FASTA format
    alignment0 = AlignIO.read(alignment_file, "fasta")
    alignment = alignment0[:, 29:-29]

    # Generate a file name for low frequency haplotypes based on the output file name
    low_seq_output_file = output_file.split(".")[0] + "_lowFreqHap.fasta"

    # Store the first sequence for reference and remove it from the analysis
    first_seq = alignment[0]
    alignment = alignment[1:]

    # Count the occurrence of each sequence in the alignment
    seq_dict = Counter(str(record.seq) for record in alignment)

    # Set the threshold for considering a sequence as a frequent haplotype
    threshold = 30
    #threshold = len(alignment) * 0.01

    # Set the lower threshold to identify low frequency haplotypes
    low_seq_threshold = 2

    # Initialize alignments for storing high and low frequency haplotypes
    new_alignment = MultipleSeqAlignment([])
    low_seq_alignment = MultipleSeqAlignment([])

    # Initialize counters for haplotypes
    hap = 0
    low_seq_hap = 0

    # Iterate through each sequence and its count
    for seq, count in seq_dict.items():
        if count >= threshold:
            # If the count is above the threshold, add to the high frequency alignment
            hap += 1
            new_record = SeqRecord(
                seq=Seq(seq), id=f"{first_seq.id}Hap-{hap}_({count})", description=""
            )
            new_alignment.append(new_record)
        elif low_seq_threshold < count < threshold:
            # If the count is between the low and high threshold, add to the low frequency alignment
            low_seq_hap += 1
            low_seq_record = SeqRecord(
                seq=Seq(seq),
                id=f"{first_seq.id}LowSeqHap-{low_seq_hap}_({count})",
                description="",
            )
            low_seq_alignment.append(low_seq_record)

    # Write the high frequency haplotypes to the output file
    AlignIO.write(new_alignment, output_file, "fasta")

    # If there are any low frequency haplotypes, write them to a separate file
    if low_seq_hap > 0:
        AlignIO.write(low_seq_alignment, low_seq_output_file, "fasta")
