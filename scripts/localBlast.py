import subprocess
import os
from Bio import SeqIO

def check_blast_db_exists(db_path, db_type):
    """
    Checks if the BLAST database files already exist.

    :param db_path: The path and base name of the BLAST database.
    :param db_type: Type of database ('nucl' for nucleotides or 'prot' for proteins).
    :return: True if the database exists, False otherwise.
    """
    # Extensions for BLAST database files
    extensions = ['.nin', '.nhr', '.nsq'] if db_type == 'nucl' else ['.pin', '.phr', '.psq']

    # Check for the existence of database files
    return all(os.path.isfile(db_path + ext) for ext in extensions)



def make_blast_db(input_file, db_type, out_db):
    """
    Creates a BLAST database using makeblastdb.

    :param input_file: Path to the input file (e.g., a FASTA file).
    :param db_type: Type of database ('nucl' for nucleotides or 'prot' for proteins).
    :param out_db: Path and name of the resulting database.
    :return: None
    """

    # Check if the BLAST database already exists
    if check_blast_db_exists(out_db, db_type):
        print(f"BLAST database {out_db} already exists. Skipping database creation.")
        return
    
    # makeblastdb command
    command = [
        'makeblastdb',
        '-in', input_file,
        '-dbtype', db_type,
        '-out', out_db
    ]

    # Execute the command
    process = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # Check for errors
    if process.returncode != 0:
        print(f"Error creating BLAST database: {process.stderr}")
    else:
        print(f"BLAST database successfully created at {out_db}")


def find_query_files(directory_path):
    """
    Find all files in the given directory that end with '_cluster.fasta'.

    :param directory_path: Path to the directory to search.
    :return: List of file paths.
    """
    query_files = []
    for file in os.listdir(directory_path):
        if file.endswith("_cluster.fasta") or file.endswith("_lowFreqHap.fasta"):
            query_files.append(os.path.join(directory_path, file))
    return query_files



def run_blast(directory_path, database, blast_type='blastn'):
    """
    Runs BLAST analysis for each file ending with '_cluster.fasta' in the given directory.
    The output files will have the same name as the query files but with '_blast.txt' extension.

    :param directory_path: Path to the directory containing query files.
    :param database: Path to the BLAST database.
    :param blast_type: Type of BLAST to run (e.g., 'blastn', 'blastp', etc.).
    :return: None
    """
    query_files = find_query_files(directory_path)

    for query_file in query_files:
        if query_file.endswith("_cluster.fasta"):
            output_file = query_file.replace("_cluster.fasta", "_cluster_blast.tsv")
            # Count the number of sequences in the query file
            sequence_count = sum(1 for _ in SeqIO.parse(query_file, "fasta"))

            if sequence_count == 1:
                # If only one sequence, run BLAST normally
                _run_single_blast(query_file, database, output_file, blast_type)
            else:
                # If more than one sequence, run BLAST for each and combine results
                combined_results = ""
                for record in SeqIO.parse(query_file, "fasta"):
                    temp_query_file = "temp_" + record.id + ".fasta"
                    temp_output_file = "temp_" + record.id + "_blast_results.txt"

                    # Write the individual sequence to a temporary file
                    SeqIO.write(record, temp_query_file, "fasta")

                    # Run BLAST for the individual sequence
                    _run_single_blast(temp_query_file, database, temp_output_file, blast_type)

                    # Read and append the result
                    with open(temp_output_file, "r") as file:
                        combined_results += file.read() + "\n\n\n"  # Add two empty lines

                    # Clean up temporary files
                    os.remove(temp_query_file)
                    os.remove(temp_output_file)

                # Write combined results to the final output file
                with open(output_file, "w") as file:
                    file.write(combined_results)
        elif query_file.endswith("_lowFreqHap.fasta"):
            output_file = query_file.replace("_lowFreqHap.fasta", "_lowFreqHap_blast.tsv")
            # Count the number of sequences in the query file
            sequence_count = sum(1 for _ in SeqIO.parse(query_file, "fasta"))

            if sequence_count == 1:
                # If only one sequence, run BLAST normally
                _run_single_blast(query_file, database, output_file, blast_type)
            else:
                # If more than one sequence, run BLAST for each and combine results
                combined_results = ""
                for record in SeqIO.parse(query_file, "fasta"):
                    temp_query_file = "temp_" + record.id + ".fasta"
                    temp_output_file = "temp_" + record.id + "_blast_results.txt"

                    # Write the individual sequence to a temporary file
                    SeqIO.write(record, temp_query_file, "fasta")

                    # Run BLAST for the individual sequence
                    _run_single_blast(temp_query_file, database, temp_output_file, blast_type)

                    # Read and append the result
                    with open(temp_output_file, "r") as file:
                        combined_results += file.read() + "\n\n\n"  # Add two empty lines

                    # Clean up temporary files
                    os.remove(temp_query_file)
                    os.remove(temp_output_file)

                # Write combined results to the final output file
                with open(output_file, "w") as file:
                    file.write(combined_results)


def _run_single_blast(query_file, database, output_file, blast_type):
    """
    Helper function to run a single BLAST search.

    :param query_file: Path to the query file.
    :param database: Path to the BLAST database.
    :param output_file: Path to the file where results will be saved.
    :param blast_type: Type of BLAST to run.
    :return: None
    """
    command = [
        blast_type,
        '-query', query_file,
        '-db', database,
        '-out', output_file,
        '-outfmt', "6 qseqid sseqid length nident pident mismatch gapopen qstart qend sstart send evalue bitscore",  # Incluyendo nident para identidades
        '-max_target_seqs', '10'
    ]

    process = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    if process.returncode != 0:
        print(f"Error running BLAST: {process.stderr}")
    else:
        # Append column names to the output
        _add_column_names_to_output(output_file)

def _add_column_names_to_output(output_file):
    """
    Helper function to add column names to the BLAST output.

    :param output_file: Path to the BLAST output file.
    :return: None
    """
    column_names = "query_id\tsubject_id\talignment_length\tn_identities\tperc_identity\tmismatches\tgap_opens\tq_start\tq_end\ts_start\ts_end\te_value\tbit_score\n"
    with open(output_file, 'r') as file:
        blast_output = file.read()
    with open(output_file, 'w') as file:
        file.write(column_names + blast_output)

# # # Example usage for make_blast_db
# make_blast_db('../blast/HmtG_database_PacBio.fasta', 'nucl', '../blast/HmtG_database_PacBio')

# # # # Example usage for run_blast
# run_blast('../AE40', '../blast/HmtG_database_PacBio')
