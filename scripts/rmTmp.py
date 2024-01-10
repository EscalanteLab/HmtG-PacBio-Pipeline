import os

def remove_temp_files(rmFiles):
    # Removing temporal files
    if rmFiles == "yes":
        for file in os.listdir():
            file_path = file
            if file.endswith("bin.txt") \
                or file.endswith("_mafftFinal.fasta")\
                or file.endswith("_corrected.fasta")\
                or file.endswith(".fa")\
                or file.endswith("-.txt")\
                or file.endswith("_mafftRaw.fasta")\
                or file.endswith("consensus.fasta-_HaplotypesFinal.fasta")\
                or file.endswith("consensus.fasta-_corrHapng.fasta")\
                or file.endswith("_corrHapng.fasta")\
                or file.endswith("_corrHap.fasta")\
                or file.endswith("_consensus.fasta-_corrHap.fasta")\
                or file.endswith("_aliHap.fasta")\
                or file.endswith("_consensus.fasta-_aliHap.fasta")\
                or file.endswith("_RawHapng.fasta")\
                or file.endswith("_consensus.fasta-_RawHapng.fasta")\
                or file.endswith("_consensus.fasta-_RawHap.fasta")\
                or file.endswith("_RawHap.fasta")\
                or file.endswith("_consensus.fasta")\
                or file.endswith("_nolong.fasta"):
                if os.path.exists(file_path):
                    os.remove(file_path)
                    print(f"Temporal file has been deleted {file_path}.")
                else:
                    print(f"Temporal file does not exist {file_path}.")
        print("Done!")
        print("HmtG pipeline has completed the analysis. Thank you!")
    else:
        print("Done!")
        print("HmtG pipeline has completed the analysis. Thank you!")
