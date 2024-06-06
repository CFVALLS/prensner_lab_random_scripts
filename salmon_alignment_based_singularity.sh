#!/bin/bash

# Define paths
SAMPLES="/data/cfvall/input/riboseq_star/sorted/bam_file_paths.txt"
CONTAINER_PATH="/data/cfvall/containers/salmon_1.5.0.sif"
REFERENCE_FASTA="/data/shared/assemblies/GENCODE_Phase1_Riboseq_ORFs/Ribo-seq_ORFs.fa"
OUTPUT_DIR_BASE="/data/cfvall/input/riboseq_star/sorted"

# Read each line from the SAMPLES file
while IFS= read -r INPUT_BAM; do
    # Extract the sample name from the BAM file path
    SAMPLE_NAME=$(basename "$INPUT_BAM" .genome.sorted.bam)
    OUTPUT_DIR="${OUTPUT_DIR_BASE}/${SAMPLE_NAME}_salmon_quant"

    echo "Processing sample: $SAMPLE_NAME"
    echo "Input BAM: $INPUT_BAM"
    echo "Output directory: $OUTPUT_DIR"

    # Run salmon quant with Singularity
    singularity exec -B /data:/data "$CONTAINER_PATH" \
        salmon quant \
        --libType A \
        -t "$REFERENCE_FASTA" \
        -a "$INPUT_BAM" \
        -o "$OUTPUT_DIR" \
        | sed "s/^/${SAMPLE_NAME} /"

done < "$SAMPLES"

echo "All samples processed."
