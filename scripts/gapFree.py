from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def remove_gap_columns(input_alignment_file, output_alignment_file, trim_at_position):
    # Parse the input alignment file
    alignments = list(SeqIO.parse(input_alignment_file, "fasta"))

    # Get the number of sequences and the alignment length
    num_sequences = len(alignments)
    alignment_length = len(alignments[0].seq)

    # Set to store the positions with gaps
    gap_positions = set()

    # Identify the positions with gaps in at least 5% of the sequences
    for i in range(alignment_length):
        gap_count = 0
        for j in range(num_sequences):
            if alignments[j].seq[i] == "-":
                gap_count += 1

        # Check if the gap count exceeds 10% of the total number of sequences
        if gap_count >= num_sequences * 0.1:
            gap_positions.add(i)

    # Create a new alignment without the gap columns
    new_alignments = []
    for alignment in alignments:
        # Generate the sequence without the gap positions
        seq = "".join(
            [alignment.seq[i] for i in range(alignment_length) if i not in gap_positions]
        )
        if trim_at_position == "yes":
            # Trim the sequence at position 5866 (adjust as needed)
            trimmed_seq = seq[:5866]

            # Create a SeqRecord for the new alignment sequence
            new_alignment = SeqRecord(Seq(trimmed_seq), id=alignment.id, description="")

            # Add the new alignment sequence to the list
            new_alignments.append(new_alignment)
        else:
            new_alignment = SeqRecord(Seq(seq), id=alignment.id, description="")
            new_alignments.append(new_alignment)

    # Write the new alignment to the output file in FASTA format
    SeqIO.write(new_alignments, output_alignment_file, "fasta")
