#!/bin/bash

# Define paths
SAMPLES="/data/cfvall/input/trimgalore/fastq_file_paths.txt"
CONTAINER_PATH="/data/cfvall/containers/salmon_1.5.0.sif"
OUTPUT_DIR_BASE="/data/cfvall/input/riboseq_salmon"
TRANSCRIPTS_INDEX="/data/shared/assemblies/GENCODE_Phase1_Riboseq_ORFs/salmon_index"

# Read each line from the SAMPLES file
while IFS= read -r INPUT_FASTQ; do
    # Extract the sample name from the FASTQ file path
    SAMPLE_NAME=$(basename "$INPUT_FASTQ" _trimmed.fq.gz)
    OUTPUT_DIR="${OUTPUT_DIR_BASE}/${SAMPLE_NAME}_salmon_quant"

    echo "Processing sample: $SAMPLE_NAME"
    echo "Input FASTQ: $INPUT_FASTQ"
    echo "Output directory: $OUTPUT_DIR"

    # Create output directory if it doesn't exist
    mkdir -p "$OUTPUT_DIR"

    # Run salmon quant with Singularity
    singularity exec -B /data:/data "$CONTAINER_PATH" \
        salmon quant \
        -i "$TRANSCRIPTS_INDEX" \
        --libType A \
        -r "$INPUT_FASTQ" \
        --validateMappings \
        -o "$OUTPUT_DIR" \
        | sed "s/^/${SAMPLE_NAME} /"
    
done < "$SAMPLES"

echo "All samples processed."
