import pandas as pd
import argparse
import os
import csv


def process_bed_line(input_line: list, df_in: pd.DataFrame) -> None:
    """
    Processes a single line from a BED file and adds it to the output DataFrame.

    Args:
        input_line: A list of strings representing the BED file line.
        df_in: The pandas DataFrame to store the processed GTF data.
    """
    try:
        block_nums = int(input_line[9])
        chrom = input_line[0]
        source = 'Gencode_orfeome_phase_1'
        score = input_line[4]
        strand = input_line[5]
        frame = '.'
        transcripts = input_line[22]
        if len(transcripts) > 20:
            transcripts = transcripts.split(',')[0]
        attributes = 'gene_id "{}|{}"; transcript_id "{}"; orf_id "{}"; gene_name "{}"; gene_type "{}"'.format(input_line[17], input_line[3], transcripts, input_line[3], input_line[18], input_line[16])
        if block_nums <= 1:
            start = int(input_line[1]) + 1  # BED is 0-based; GTF is 1-based
            end = int(input_line[2])
            attributes = attributes + '; exon_number 1'
            new_line = [chrom, source, "CDS", start, end, score, strand, frame, attributes]
            df_in.loc[len(df_in)] = new_line
        else:
            block_sizes = [int(x) for x in input_line[10].split(',') if x]
            block_starts = [int(x) for x in input_line[11].split(',') if x]

            block_iter_list = list(range(1,block_nums+1))
            if strand == '-':
                block_iter_list = block_iter_list[::-1]

            loop_count = 0
            for size, start in zip(block_sizes, block_starts):
                exon_id = block_iter_list[loop_count]
                new_start = int(input_line[1]) + start + 1  # BED is 0-based; GTF is 1-based
                new_end = new_start + size - 1
                # attributes = attributes + '; exon_id {}'.format(exon_id)
                new_line = [chrom, source, "CDS", new_start, new_end, score, strand, frame, attributes + '; exon_id {}'.format(exon_id)]
                df_in.loc[len(df_in)] = new_line
                loop_count += 1
    
    except ValueError as e:
        print(f"Error processing line (invalid integer format): {input_line}")
        print(f"Error details: {e}")
    except Exception as e:
        print(f"Error processing line: {input_line}")
        print(f"Error details: {e}")

def main(bed_file, output_dir):
    try:
        # Initialize an empty DataFrame for the output GTF file
        gtf_columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
        gtf_df = pd.DataFrame(columns=gtf_columns)

        # Load the BED file
        gencode_orfs = pd.read_csv(bed_file, sep=' ', header=None)

        # Process each line in the BED file
        gencode_orfs.apply(lambda row: process_bed_line(row, gtf_df), axis=1)

        # Extract the base filename (without extension)
        filename, _ = os.path.splitext(os.path.basename(bed_file))
        output_file = filename + '.gtf'

        if output_dir:  # If an output directory is provided
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            output_file = os.path.join(output_dir, output_file)

        # Save the resulting GTF DataFrame to a file
        gtf_df.to_csv(output_file, sep='\t', header=False, index=False, quoting=csv.QUOTE_NONE)
    except Exception as e:
        print("Failed to process the file: {}".format(e))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert BED file to GTF format.')
    parser.add_argument('bed_file', type=str, help='Path to the BED file')
    parser.add_argument('--output_dir', type=str, default='gtf', help='Directory to save the output GTF file')

    args = parser.parse_args()

    if len(args.output_dir) == 0:
        main(args.bed_file)
    else:
        main(args.bed_file, args.output_dir)

