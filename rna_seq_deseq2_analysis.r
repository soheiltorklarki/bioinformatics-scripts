# Title: Differential Gene Expression Analysis Pipeline
# Author: Soheil Tork Larki
# Description: Standard DESeq2 workflow for RNA-Seq data processing in R

library(DESeq2)
library(ggplot2)

# 1. Simulate Expression Count Data
set.seed(123)
cts <- matrix(rpois(1000, lambda = 100), nrow = 100, ncol = 10)
rownames(cts) <- paste0("Gene_", 1:100)
colnames(cts) <- paste0("Sample_", 1:10)

# 2. Define Sample Metadata
col_data <- data.frame(
  condition = factor(rep(c("Control", "Treated"), each = 5)),
  row.names = colnames(cts)
)

# 3. Create DESeq2 Dataset
dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = col_data,
  design = ~ condition
)

# 4. Run Differential Expression Analysis
dds <- DESeq(dds)
res <- results(dds)

# 5. Filter Significant Genes (p-adj < 0.05)
sig_genes <- subset(res, padj < 0.05)

# Output summary
cat("Pipeline finished. Number of differentially expressed genes:", nrow(sig_genes), "\n")