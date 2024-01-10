from Bio import SeqIO
import os
import itertools
import pandas as pd
import subprocess

def p_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("The sequences must have the same length in order to calculate the p-distance.")

    diff = 0
    total = 0

    for a, b in zip(seq1, seq2):
        if a != '-' and b != '-' and a != 'N' and b != 'N':
            total += 1
            if a != b:
                diff += 1

    return diff / total if total > 0 else None

def run_mafft(input_file_path, output_file_path):
    with open(output_file_path, 'w') as outfile:
        subprocess.run(["mafft", input_file_path], stdout=outfile)

def calculate_pairwise_distance(file_path):
    records = list(SeqIO.parse(file_path, "fasta"))
    matrix = pd.DataFrame(index=[rec.id for rec in records], columns=[rec.id for rec in records])

    for i, rec1 in enumerate(records):
        for j, rec2 in enumerate(records):
            if i <= j:
                dist = p_distance(str(rec1.seq), str(rec2.seq))
                matrix.at[rec1.id, rec2.id] = dist
                matrix.at[rec2.id, rec1.id] = dist
    return matrix

def calculate_between_files_distance(file1, file2):
    records1 = list(SeqIO.parse(file1, "fasta"))
    records2 = list(SeqIO.parse(file2, "fasta"))

    with open("temp_sequences.fasta", "w") as f:
        SeqIO.write(records1 + records2, f, "fasta")

    run_mafft("temp_sequences.fasta", "temp_sequences_aligned.fasta")
    records = list(SeqIO.parse("temp_sequences_aligned.fasta", "fasta"))

    aligned_records1 = records[:len(records1)]
    aligned_records2 = records[len(records1):]

    matrix = pd.DataFrame(index=[rec.id for rec in aligned_records1], columns=[rec.id for rec in aligned_records2])

    for rec1 in aligned_records1:
        for rec2 in aligned_records2:
            dist = p_distance(str(rec1.seq), str(rec2.seq))
            matrix.at[rec1.id, rec2.id] = dist

    os.remove("temp_sequences.fasta")
    os.remove("temp_sequences_aligned.fasta")

    return matrix


def calculate_all_distances():
    fasta_files = [file for file in os.listdir('.') if file.endswith('_cluster.fasta')]
    output_str = ""

    # Within-file distances
    for file in fasta_files:
        matrix = calculate_pairwise_distance(file)
        output_str += f"Within-cluster p-distance for {file}\n"
        output_str += "="*40 + "\n"
        output_str += matrix.to_string() + "\n\n"

    # Between-cluster distances
    for file1, file2 in itertools.combinations(fasta_files, 2):
        matrix = calculate_between_files_distance(file1, file2)
        output_str += f"Between-cluster p-distance for {file1} and {file2}\n"
        output_str += "="*40 + "\n"
        output_str += matrix.to_string() + "\n\n"

    return output_str
