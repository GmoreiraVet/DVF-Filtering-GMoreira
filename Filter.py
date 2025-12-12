import os
import pandas as pd
import subprocess
from Bio import SeqIO

# ==========================
# PARA EDITAR LOCALIZAÇÕES DOS FICHEIROS
# ==========================
fastq_file = "/path/to/your/ont_reads.fastq"   # single-end ONT reads
output_folder = "/path/to/output_folder"
dvf_path = "/path/to/dvf.py"                  # path to dvf.py in your conda env
cutoff_length = 500                            # minimum contig/read length
core_num = 4                                   # number of CPU cores to use

strict_score = 0.9
lenient_score = 0.75
pval_threshold = 0.05

# ==========================
# CRIAR PASTA DE OUTPUT
# ==========================
os.makedirs(output_folder, exist_ok=True)

# ==========================
# FASTQ > FASTA
# ==========================
fasta_file = os.path.join(output_folder, "ont_reads.fasta")

with open(fasta_file, "w") as fasta_out:
    for record in SeqIO.parse(fastq_file, "fastq"):
        SeqIO.write(record, fasta_out, "fasta")

# ==========================
# CORRER DEEPVIRFINDER
# ==========================
subprocess.run([
    "python", dvf_path,
    "-i", fasta_file,
    "-o", output_folder,
    "-l", str(cutoff_length),
    "-c", str(core_num)
])

# DVF output
dvf_output_file = os.path.join(output_folder, "ont_reads.fasta.score")

# ==========================
# LER OS RESULTADOS
# ==========================
df = pd.read_csv(dvf_output_file, sep="\t")

# ==========================
# FILTRAR COM OS DOIS GRUPOS DE PARAMETROS
# LENIENT E STRICT
# ==========================
strict_contigs = set(df.loc[(df['score'] >= strict_score) & (df['pvalue'] <= pval_threshold), 'contig'])
lenient_contigs = set(df.loc[(df['score'] >= lenient_score) & (df['pvalue'] <= pval_threshold), 'contig'])

# ======================================================
# CRIAR NOVOS FASTQ COM OS FILTROS APLICADOS
# ======================================================
def filter_fastq(input_fastq, output_fastq, contig_set):
    with open(output_fastq, 'w') as out:
        for record in SeqIO.parse(input_fastq, "fastq"):
            rid = record.id.split()[0]
            if rid in contig_set:
                SeqIO.write(record, out, "fastq")

strict_fastq = os.path.join(output_folder, "viral_reads_strict.fastq")
lenient_fastq = os.path.join(output_folder, "viral_reads_lenient.fastq")

filter_fastq(fastq_file, strict_fastq, strict_contigs)
filter_fastq(fastq_file, lenient_fastq, lenient_contigs)

print("Done! Filtered FASTQ files saved in:", output_folder)
