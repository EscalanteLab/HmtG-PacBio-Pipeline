import os
import subprocess


def concat(filesC0, filesC1):
    for i in filesC0:
        for j in filesC1:
            if i.split("-")[0] == j.split("-")[0]:
                outf = i.split("-")[0] + "-_RawHap.fasta"
                with open(outf, "w") as f:
                    cat_command = ["cat", i, j]
                    subprocess.run(cat_command, stdout=f)
