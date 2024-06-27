cran_mirror <- "https://cran.rstudio.com/"

# Install and load necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = cran_mirror)
}

# Ensure BiocManager is loaded
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    stop("BiocManager not available.")
}


BiocManager::install("SummarizedExperiment")
BiocManager::install("DESeq2")

library(readr)
library(dplyr)
library(DESeq2)


# OPTION_1
# import and parse files
count_path <- '/data/cfvall/input/RNAseq_DE/Nextflow_analysis/salmon.merged.gene_counts.tsv'
count <- read_delim(count_path, delim='\t', escape_double=F, trim_ws=T)
# rename ugly column names 
colnames(count) <- sub("//_RNA*", "", colnames(count)) 
# OPTION_2
# import RDS file to create a SummerizedExperiment object 
rds_path <- '/data/cfvall/input/RNAseq_DE/Nextflow_analysis/salmon.merged.gene_counts.rds'
count.SE <- readRDS(rds_path)

# # Create a DDS object using S.E object
# ddsSE <- DESeqDataSet(se, design = ~ cell + dex)

# # Filter-Out low count rows
# smallestGroupSize <- 3
# keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
# dds <- dds[keep,]



# Assign factors to each sample 
# ew vs no_ew
# ew: CADOES1, A673, RDES 
# no_we: THP1, A375, A549, Jurkat, HEPG2
colData = 




## Create DSeq2 object



