#!/bin/bash

# >>>>>> this one is working: 
# bedtools intersect -a /home/cfvall/data/projects/2024_variants_ncORF/gnomad.wgs.y.vcf -b /home/cfvall/shared/assemblies/GENCODE_Phase1_Riboseq_ORFs/Ribo-seq_ORFs.bed -header > testclear
# Note that it would exclude everything outside that regions, maybe its better to construct an annotation file with the parental genes


# Use bedtools to intersect VCF file with BED file and obtain a reduced VCF with only overlapping sequences.
vcf_file="/data/cfvall/projects/2024_variants_ncORF/gnomad.wgs.y.vcf"
bed_file="/home/cfvall/shared/assemblies/GENCODE_Phase1_Riboseq_ORFs/Ribo-seq_ORFs_nochr.bed.gz"
output_file="${vcf_file%.vcf}.intersected.ncORFs.vcf"

# # Ensure the BED file is indexed
# if [ ! -f "${bed_file}.gz" ]; then
#     echo "Bgzipping the BED file..."
#     bgzip -c "$bed_file" > "${bed_file}.gz"
# fi

# if [ ! -f "${bed_file}.gz.tbi" ]; then
#     echo "Indexing the gzipped BED file..."
#     tabix -p bed "${bed_file}.gz"
# fi

# Intersect the VCF file with the BED file
echo "Intersecting VCF file with BED file..."
bedtools intersect -a "$vcf_file" -b "${bed_file}.gz" -header > "$output_file"

echo "Intersection complete. Output written to $output_file"