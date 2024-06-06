import os
import subprocess
import pandas as pd

input_fasta = "/home/cfvall/input/query.fa"
database_path = os.path.expanduser("/home/cfvall/shared/assemblies/GRCh38_MANE_transcripts/MANEdb")
output_file = '/home/cfvall/output/result_2.tsv'

blast_command = [
    'blastn',
    '-query', input_fasta,
    '-task', 'blastn',
    '-db', database_path,
    '-evalue', '10',
    '-word_size', '4',
    '-num_threads', '14',
    '-out', output_file,
    '-max_target_seqs', '5',
    '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand'
]

subprocess.run(blast_command)

def read_blast_output(output_file):
    """Reads BLAST output (outfmt 6) including the strand and returns a pandas dataframe."""
    return pd.read_csv(output_file, sep="\t",
                       names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sstrand"],  # Include sstrand
                       index_col="qseqid")

# Read and print the BLAST output dataframe
df = read_blast_output(output_file)
print(df.head())
