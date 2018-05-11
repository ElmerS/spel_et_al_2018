#######################################################################
# DESeq2 workflow for analysis of CRISPR Screen in Spel et al. (2018) #
# Elmer Stickel, Netherlands Cancer Institute, 2018.                  #
# Contact: elmer.stickel [at] posteo.net                              #
# Based on: https://www.bioconductor.org/help/workflows/rnaseqGene/   #
#######################################################################


# Load packages
source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2") # Install DESeq2 if it wasn't been installed before
library("DESeq2")
#biocLite("IRanges") # Install IRanges if it wasn't been installed before
library("IRanges")

# Read the data
guides_table <- read.table('unique_guides_underscore.txt', col.names='sgRNA')

pos_rep_1_table <- read.table("sample_w.counts", header=F, col.names = c('P1', 'sgRNA')) # Counts of sorted population with high NFkB and high MHC, replicate 1
pos_rep_2_table <- read.table("sample_x.counts", header=F, col.names = c('P2', 'sgRNA')) # Counts of sorted population with high NFkB and high MHC, replicate 2
neg_rep_1_table <- read.table("sample_y.counts", header=F, col.names = c('N1', 'sgRNA')) # Counts of sorted population with low NFkB and low MHC, replicate 1
neg_rep_2_table <- read.table("sample_z.counts", header=F, col.names = c('N2', 'sgRNA')) # Counts of sorted population with low NFkB and low MHC, replicate 2

# Merge the counts from all four populations into a single dataframe
countdata_df <- Reduce(function(...) merge(..., all=TRUE, by='sgRNA'), list(guides_table, pos_rep_1_table, pos_rep_2_table, neg_rep_1_table, neg_rep_2_table))
# Change to rownames to sgRNAs
rownames(countdata_df) <- countdata_df$sgRNA

# Set rows with NAs to zero
countdata_df[is.na(countdata_df)] <- 0

countdata <- data.matrix(within(countdata_df, rm('sgRNA')))
population <- factor(c(rep("pos", 2), rep("neg", 2)))
coldata <- data.frame(row.names=colnames(countdata), population)
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~population)
dds <- DESeq(dds)
# dds <- dds[ rowSums(counts(dds)) > 20, ] # Optional: use to remove guides where the sum of counts across all populations is 20 or less.
res <- results(dds) # Get the results

# These are the columns:
mcols(res, use.names = TRUE)

# Apply Log Fold Change shrinkage
res <- lfcShrink(dds, contrast=c("population","pos","neg"), res=res)

# Order values on FDR-adjusted p-values
resOrdered <- res[order(res$padj),]

# Make a dataframe of the ordered data
resOrderedDF <- as.data.frame(resOrdered)

# Export results to results.csv
write.csv(resOrderedDF, file = "intermediate_results.csv")

print("Completed! Continue with postprocessing")
