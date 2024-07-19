#!/bin/bash

vcf_file="/data/cfvall/projects/2024_variants_ncORF/gnomad.wgs.y.vcf"
bed_file="/data/shared/assemblies/GENCODE_Phase1_Riboseq_ORFs/Ribo-seq_ORFs_nochr.bed.gz"
gtf_file="/data/shared/assemblies/GENCODE_Phase1_Riboseq_ORFs/ncORFs_complete.gtf.gz"
genome_fasta="/data/shared/assemblies/GRCh38_110/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
fasta="/data/shared/assemblies/GRCh38_vep/homo_sapiens/112_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
vep_cache_dir="/data/shared/assemblies/GRCh38_vep"

# Ensure the reference genome FASTA file is indexed
if [ ! -f "${genome_fasta}.fai" ]; then
    echo "Indexing the reference genome FASTA file..."
    samtools faidx "$genome_fasta"
fi

# Ensure the BED file is indexed
if [ ! -f "${bed_file}.tbi" ]; then
    echo "Indexing the BED file..."
    tabix -p bed "$bed_file"
fi

# Create output directory if it doesn't exist
mkdir -p vep_output

# Run VEP in offline mode with multiple forks
vep \
    -i "$vcf_file" \
    -v \
    --offline \
    --dir_cache "$vep_cache_dir" \
    --custom "$gtf_file,Ribo-seq_ORFs,gtf,overlap" \
    --fasta "$fasta" \
    --fork 10 \
    --output_file vep_output/annotated_variants.txt
