from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

def correct_sequences(input_msa_file, output_msa_file):
    alignment = AlignIO.read(input_msa_file, "fasta")
    corrected_alignment = []

    # Consensus sequence is assumed to be the first sequence in the alignment
    consensus_seq = alignment[0]
    corrected_alignment.append(consensus_seq)

    for record in alignment[1:]:
        corrected_seq_str = ''

        for cons, seq in zip(str(consensus_seq.seq).upper(), str(record.seq).upper()):
            if cons != 'N':
                corrected_seq_str += cons
            else:
                corrected_seq_str += seq

        corrected_seq = SeqRecord(Seq(corrected_seq_str), id=record.id)
        corrected_alignment.append(corrected_seq)

    corrected_alignment = MultipleSeqAlignment(corrected_alignment)
    AlignIO.write(corrected_alignment, output_msa_file, "fasta")
