import sys
import os
import argparse
import subprocess

sys.path.append('scripts/')

## import functions
from readsRedirection import primerDetection
from mafft import mafftRaw, mafftFinal, mafftHap
from fasta2binary import fasta2bin
from VAEmethod import VAE_model
from DBScan import extract_cluster_labels_dbscan
from cluster2fasta import txt2fasta
from consensus import cons
from rmTmp import remove_temp_files
from concatenate import concat
from corrHap import correct_sequences
from removeGap import remove_gaps_from_multifasta
from gapFree import remove_gap_columns
from countHap import countHaplotypes
from localBlast import make_blast_db, run_blast
import summary


def main():
    ## Arguments 
    parser = argparse.ArgumentParser(description='Haemosporidian Mitochondrial Genome PacBio Pipeline (HmtG-PacBio Pipeline)')
    ### main arguments
    parser.add_argument('-rR', '--rawReads', type=str, help="rawReads from PacBio sequencing", required=True)
    parser.add_argument('-pF', '--primerF', type=str, help="primerF 5'-3', default=GATTCTCTCCACACTTCAATTCGTACTTC", default='GATTCTCTC')
    parser.add_argument('-pR', '--primerR', type=str, help="primerR 3'-5', default=GAAGTACGAATTGAAGTGTGGAGAGAATC", default='AGACCGAACCTTGGACTC')
    parser.add_argument('-eps', '--epsDBScan', type=float, help="epsilon value for DBScan clustering, default=1.0", default=1.0)
    parser.add_argument('-rF', '--RemoveFiles', type=str, help="Removing temporal files, default=yes", default='yes')
    parser.add_argument('-rB', '--blastn', type=str, help="Run blastn locally, default=yes", default='yes')

    args = parser.parse_args()

    ## files
    rawReads = args.rawReads
    nameSample = str(rawReads).split(".")[0]
    dirSample = str(rawReads).split("/")[0]
    dirb = os.getcwd()
    dirc = os.getcwd()+"/"+dirSample

    ## run readsRedirection
    primerF = args.primerF
    primerR = args.primerR
    primerDetection(rawReads, primerF, primerR)

    rmFiles = args.RemoveFiles
    
    ## run raw aligment

    no_long_exists = any(f.endswith("_mafftRaw.fasta") for f in os.listdir(dirc))
    if not no_long_exists:
        print("Running first alignment...")
        mafftRawOut = str(rawReads).split(".")[0] + "_nolong.fasta"
        mafftRaw(mafftRawOut, dirc)
    else:
        mafftRawOut = str(rawReads).split(".")[0] + "_mafftRaw.fasta"
        print("A file with extension '_mafftRaw.fasta' was found, no initial alignment is required.")

    ## run convert DNA seq to Binary format

    no_long_exists = any(f.endswith("_bin.txt") for f in os.listdir(dirc))
    if not no_long_exists:
        fasta2binOut= str(rawReads).split(".")[0] + "_mafftRaw.fasta"
        fasta2bin(fasta2binOut, dirc)
        print("Converting DNA into bytes...")
    else:
        print("A file with extension '_bin.txt' was found, no convertion is required.")

    ## run VAE program and clustering
    epsilon = args.epsDBScan
    VAErunOut = str(rawReads).split(".")[0] + "_bin.txt"
    mu, VAErunData = VAE_model(VAErunOut, epsilon)
    nameCluster = str(rawReads).split(".")[0]
    extract_cluster_labels_dbscan(mu, VAErunData, nameCluster, epsilon)

    ## run cluster2fasta
    os.chdir(dirc)
    txtCluster = [f for f in os.listdir() if f.endswith('-.txt') and f.startswith(nameSample.split("/")[0])]
    for headers_file in txtCluster:
        mapped_output_file = headers_file.split(".")[0] + ".fa"
        txt2fasta(mafftRawOut.split("/")[1], headers_file, mapped_output_file)

    print("Running last alignment...")
    ## run final aligment
    filesM = [f for f in os.listdir() if f.endswith('.fa') and f.startswith(nameSample.split("/")[0])]
    mafftFinal(filesM)

    ## consensus output
    filesC = [f for f in os.listdir(dirc) if f.endswith('_mafftFinal.fasta') and f.startswith(nameSample.split("/")[0])]
    cons(filesC)

    ## concatenate output
    filesC0  = [f for f in os.listdir() if f.endswith('_mafftFinal.fasta') and f.startswith(nameSample.split("/")[0])]
    filesC1  = [f for f in os.listdir() if f.endswith('_consensus.fasta')  and f.startswith(nameSample.split("/")[0])]
    concat(filesC1, filesC0)

    ## gapfree output
    filesG0  = [f for f in os.listdir() if f.endswith('-_RawHap.fasta')]
    for file in filesG0:
        outf = str(file).split("-")[0] + "-_RawHapng.fasta"
        remove_gaps_from_multifasta(file, outf)

    ### mafft gapfree
    filesM1 = [f for f in os.listdir() if f.endswith('-_RawHapng.fasta')]
    mafftHap(filesM1)

    ## correcting Haplotype alignment
    filesHc = [f for f in os.listdir() if f.endswith('-_aliHap.fasta')]
    for file in filesHc:
        outf = str(file).split("-")[0] + "-_corrHap.fasta"
        correct_sequences(file, outf)

    ## Haplotypes gapfree output
    filesG1  = [f for f in os.listdir() if f.endswith('-_corrHap.fasta')]
    for file in filesG1:
        outf = str(file).split("-")[0] + "-_corrHapng.fasta"
        remove_gap_columns(file, outf, "no")

    ## counting Haplotypes
    filesH = [f for f in os.listdir() if f.endswith('-_corrHapng.fasta')]
    for file in filesH:
        outf = str(file).split("-")[0] + "_cluster.fasta"
        countHaplotypes(file, outf)

    ## distances summary
    dinstances = summary.calculate_all_distances()
    with open("summary_distances.tsv", "w") as f:
        f.write(dinstances)

    ## run local Blast
    if args.blastn == 'yes' or args.blastn == 'y' or args.blastn == 'Yes' or args.blastn == 'Y':
        print("Running local blastn...")
        make_blast_db(dirb + '/blast/HmtG_database_PacBio.fasta', 'nucl', dirb + '/blast/HmtG_database_PacBio')
        run_blast(dirc, dirb + '/blast/HmtG_database_PacBio')


    # removing temporal files
    remove_temp_files(rmFiles)

if __name__ == '__main__':
    main()
