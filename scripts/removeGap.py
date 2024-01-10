from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def remove_gaps_from_multifasta(input_file, output_file):
    # Read the input multifasta file
    records = SeqIO.parse(input_file, "fasta")

    # Create a list to store the gap-removed sequence records
    gap_removed_records = []

    # Iterate over each record in the multifasta file
    for record in records:
        # Get the sequence from the record and remove the gaps
        sequence = str(record.seq)  # Convert Seq to a string
        sequence = sequence.replace("-", "")

        # Create a new gap-removed record
        gap_removed_record = SeqRecord(Seq(sequence), id=record.id, description="")

        # Add the gap-removed record to the list
        gap_removed_records.append(gap_removed_record)

    # Write the gap-removed records to a new multifasta file
    with open(output_file, "w") as handle:
        SeqIO.write(gap_removed_records, handle, "fasta")
