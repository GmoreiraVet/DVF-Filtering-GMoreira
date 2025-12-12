import os
import pandas as pd
import subprocess
from Bio import SeqIO
import sys

# ==========================
# USER-EDITABLE PATHS
# ==========================
fastq_file = "/home/viroicbas2023/Documents/Gmoreira/Carracas_ticks/Carracas_Algarve/BEZ_ROG_3010/Tacheng/Dez2025/Reads/N56VG9_1_ROG-F.fastq"
output_folder = "/home/viroicbas2023/Documents/Gmoreira/Carracas_ticks/Carracas_Algarve/BEZ_ROG_3010/Tacheng/Dez2025/DVF/FILTRO"
dvf_path = "/home/viroicbas2023/Documents/DeepVirFinder/dvf.py"

cutoff_length = 500
core_num = 30

strict_score = 0.9
lenient_score = 0.75
pval_threshold = 0.05

# ==========================
# CREATE OUTPUT FOLDER
# ==========================
os.makedirs(output_folder, exist_ok=True)

# ==========================
# FASTQ → FASTA
# ==========================
fasta_file = os.path.join(output_folder, "ont_reads.fasta")
with open(fasta_file, "w") as fasta_out:
    for record in SeqIO.parse(fastq_file, "fastq"):
        SeqIO.write(record, fasta_out, "fasta")

# ==========================
# RUN DEEPVIRFINDER
# ==========================
subprocess.run([
    "python", dvf_path,
    "-i", fasta_file,
    "-o", output_folder,
    "-l", str(cutoff_length),
    "-c", str(core_num)
])

# ==========================
# DETERMINE DVF OUTPUT FILE
# ==========================
# Preferred output
dvf_output_gt500 = os.path.join(output_folder, "ont_reads.fasta_gt500bp_dvfpred.txt")
# Fallback output
dvf_output_all = os.path.join(output_folder, "ont_reads.fasta_dvfpred.txt")

if os.path.exists(dvf_output_gt500):
    dvf_output_file = dvf_output_gt500
elif os.path.exists(dvf_output_all):
    dvf_output_file = dvf_output_all
else:
    print("ERROR: DVF output was not found in:", output_folder)
    sys.exit(1)

print("Using DVF output:", dvf_output_file)

# ==========================
# READ DVF RESULTS
# ==========================
df = pd.read_csv(dvf_output_file, sep="\t")

# DVF column name is usually "name", not "contig"
contig_col = "name" if "name" in df.columns else "contig"

# ==========================
# FILTER STRICT & LENIENT
# ==========================
strict_contigs = set(df.loc[
    (df["score"] >= strict_score) &
    (df["pvalue"] <= pval_threshold),
    contig_col
])

lenient_contigs = set(df.loc[
    (df["score"] >= lenient_score) &
    (df["pvalue"] <= pval_threshold),
    contig_col
])

# ==========================
# FILTER READS BACK TO FASTQ
# ==========================
def filter_fastq(input_fastq, output_fastq, contig_set):
    with open(output_fastq, "w") as out:
        for record in SeqIO.parse(input_fastq, "fastq"):
            rid = record.id.split()[0]
            if rid in contig_set:
                SeqIO.write(record, out, "fastq")

strict_fastq = os.path.join(output_folder, "viral_reads_strict.fastq")
lenient_fastq = os.path.join(output_folder, "viral_reads_lenient.fastq")

filter_fastq(fastq_file, strict_fastq, strict_contigs)
filter_fastq(fastq_file, lenient_fastq, lenient_contigs)

print("Done! Filtered FASTQ files written to:", output_folder)
