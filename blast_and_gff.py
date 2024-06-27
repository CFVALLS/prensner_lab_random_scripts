import sys
import os
import subprocess
import pandas as pd
from Bio.Blast import NCBIXML
from Bio import SeqIO
import sys
import os
import subprocess

def run_blast(database_path, input_fasta, outputfile):
    """Run BLAST and return results in tabular format"""
    blast_command = [
        'blastn',
        '-evalue', '10',
        '-word_size', '4',
        '-num_threads', '14',
        '-db', database_path,
        '-query', input_fasta,
        '-out', outputfile,
        '-outfmt', '6',  # Tabular format
        '-max_target_seqs', '10',
    ]

    subprocess.run(blast_command, check=True)

def read_blast_output(tsv_in, gff_out, verbose=True):
    """Parse BLAST results into a GFF file"""
    with open(tsv_in) as tsv, open(gff_out, 'w') as gff:
        for line in tsv:
            fields = line.strip().split('\t')
            query_id = fields[0]
            subject_id = fields[1]
            percent_identity = fields[2]
            alignment_length = fields[3]
            mismatches = fields[4]
            gap_opens = fields[5]
            query_start = fields[6]
            query_end = fields[7]
            subject_start = fields[8]
            subject_end = fields[9]
            e_value = fields[10]
            bit_score = fields[11]
            strand = '+' if int(subject_start) < int(subject_end) else '-'
            
            gff.write(f"{subject_id}\tBLAST\tmatch\t{subject_start}\t{subject_end}\t{e_value}\t{strand}\t.\tID={query_id};Name={subject_id};E-value={e_value};Bitscore={bit_score}\n")

def test_blast_single_sequence(sequence: str, database_path, verbose=True):
    """Test function to run BLAST on a single sequence, prints the output"""
    test_fasta = "test_sequence.fasta"
    blast_output_tsv = "blast_results_test.tsv"  # Temporary file to store BLAST results
    gff_output = "test_output.gff"

    if verbose == True:
        print('Running Test')

    try:
        with open(test_fasta, 'w') as fasta_file:
            fasta_file.write(f">TestSequence\n{sequence}\n")

        run_blast(database_path=database_path, input_fasta=test_fasta, outputfile=blast_output_tsv)
        read_blast_output(tsv_in=blast_output_tsv, gff_out=gff_output)

        # Optionally, print the GFF output for inspection
        with open(gff_output, 'r') as gff_file:
            print(gff_file.read())
    except:
        print('Error in Test')
    
    # Optionally, remove the temporary files
    # os.remove(test_fasta)
    # os.remove(blast_output_tsv)
    # os.remove(gff_output)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print('Usage: python blast_to_gff.py <fasta_path> <db_path> <gff_out>')
        sys.exit(1)

    input_fasta = sys.argv[1]
    database_path = sys.argv[2]
    output_gff = sys.argv[3]

    if input_fasta == "test":
        test_sequence = "TAGGATTGTGACATATATC"  # Example test sequence
        test_blast_single_sequence(sequence=test_sequence, database_path=database_path)
    else:
        blast_output_tsv = "blast_results.tsv"  # Temporary file to store BLAST results

        run_blast(database_path=database_path, input_fasta=input_fasta, outputfile=blast_output_tsv)
        read_blast_output(tsv_in=blast_output_tsv, gff_out=output_gff)

        # Optionally, remove the temporary BLAST TSV file
        os.remove(blast_output_tsv)

