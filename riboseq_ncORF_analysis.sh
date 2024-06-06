#!/bin/bash

# Function to display usage
usage() {
    echo "Usage: $0 --input path_to_samples.csv --bed path_to_bed_file"
    exit 1
}

# Check for the correct number of arguments
if [[ $# -ne 4 ]]; then
    usage
fi

# Parse arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --input)
            samples_file="$2"
            shift # past argument
            shift # past value
            ;;
        --bed)
            bed_file="$2"
            shift # past argument
            shift # past value
            ;;
        *)
            usage
            ;;
    esac
done

# Check if the input files exist
if [[ ! -f $samples_file ]]; then
    echo "Error: Samples file '$samples_file' not found." >&2
    exit 1
fi

if [[ ! -f $bed_file ]]; then
    echo "Error: BED file '$bed_file' not found." >&2
    exit 1
fi

# Initialize an empty string for BAM files
bam_files=()

# Create an output file for combined samtools view output
combined_output="samtools_view_output.txt"
> "$combined_output" # Clear the file if it already exists

# Loop through each line in the samples CSV file
while IFS=, read -r sample_name bam_path; do
    # Skip header line
    if [[ $sample_name == "sample" ]]; then
        continue
    fi

    # Check if the BAM file exists
    if [[ ! -f $bam_path ]]; then
        echo "Error: BAM file '$bam_path' for sample '$sample_name' not found." >&2
        continue
    fi

    # Append BAM file to the list
    bam_files+=("$bam_path")

    # Run samtools commands and capture output directly into variables
    echo "Sample: $sample_name" >> "$combined_output"
    total_count=$(samtools view -c -q 30 -F 3844 "$bam_path")
    echo "Total count = $total_count" >> "$combined_output"
    in_target_count=$(samtools view -c -q 30 -F 3844 -L "$bed_file" "$bam_path")
    echo "In-target count = $in_target_count" >> "$combined_output"
    echo "" >> "$combined_output" # Add a blank line for readability

    # Check if samtools exited successfully
    if [[ $? -ne 0 ]]; then
        echo "Error: samtools command failed for sample '$sample_name'. Exiting." >&2
        exit 1
    fi

done < "$samples_file"

# Check if bedtools command is executable before attempting to run it
if ! command -v bedtools &> /dev/null; then
    echo "Error: bedtools command not found. Exiting." >&2
    exit 1
fi

# Run bedtools coverage and stream output to a file
bedtools coverage -b "${bam_files[@]}" -a "$bed_file" > multicov.txt

# Check if bedtools exited successfully
if [[ $? -ne 0 ]]; then
    echo "Error: bedtools command failed. Exiting." >&2
    exit 1
fi

echo "Combined samtools view output saved to '$combined_output'"
echo "Bedtools coverage output saved to 'multicov.txt'"

