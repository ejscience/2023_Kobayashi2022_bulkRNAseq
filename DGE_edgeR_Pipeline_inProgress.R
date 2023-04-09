#--- set working directory
setwd("/mnt/DATA/LairdLab/ernesto/kob2022_ePGCtoEGC")
#---- load packages
{
  library(dplyr)
  library(DESeq2)
  library(edgeR)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(GEOquery)
  library(tidyr)
  library(EnhancedVolcano)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(enrichplot)
}
counts_data <- read.table(file = "data/bulk/GSE174465_Normalized_counts.txt", 
                  sep = "\t", header = T, stringsAsFactors = FALSE)
sample_info <- metadata %>%
  select(1,8,44,45,46) %>%
  rename_with(.cols = 1, ~"sample_id") %>%
  rename_with(.cols = 2, ~"group") %>%
  rename_with(.cols = 3, ~"clone") %>%
  rename_with(.cols = 4, ~"cell_info") %>%
  rename_with(.cols = 5, ~"sex") %>%
  mutate(group = gsub("human ", "", group))
head(sample_info)


head(dat)
counts_data <- dat
rownames(counts_data) <- counts_data$X
counts_data$X <- NULL
head(counts_data)

if (!all(colnames(counts_data) %in% sample_info$sample_id)) {
  stop("The sample names in the counts data do not match the sample information.")
}


# Create a DGEList object from the counts data
dge <- DGEList(counts = counts_data)

# Assign group information to the DGEList object
dge$samples$group <- factor(sample_info$group)

# Estimate common dispersion
dge <- estimateCommonDisp(dge)

# Estimate tagwise dispersion
dge <- estimateTagwiseDisp(dge)

# Perform the exact test for differential expression
et <- exactTest(dge, pair=cbind("PGCLC", "EG")) # Replace "group1" and "group2" with your specific group names

# Extract and sort the results by FDR (False Discovery Rate)
degs <- topTags(et, n = Inf, p.value = 0.05, adjust.method = "BH")

# Save the differentially expressed genes to a file
write.csv(degs, file = "bulk_analysis/results/differentially_expressed_genes_edgeR.csv")



--------

# Load required libraries
library(pheatmap)

# Filter the top N most variable genes for the heatmap (adjust N as desired)
N <- 100
most_variable_genes <- head(order(rowVars(dge$counts), decreasing = TRUE), n = N)

# Extract the log-CPM values of the top N most variable genes
log_cpm <- cpm(dge, log = TRUE)
log_cpm_top_genes <- log_cpm[most_variable_genes, ]

# Create the heatmap
pheatmap(log_cpm_top_genes, 
         scale = "row", # Scale genes (rows) to have mean zero and standard deviation one
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = FALSE,
         show_colnames = TRUE)


# Load required libraries
library(ggplot2)

# Extract the results from the edgeR output
results <- et$table
results$Gene <- rownames(results)

results$FDR <- p.adjust(results$PValue, method = "BH")

# Calculate log2 fold change and -log10 adjusted p-value
results$log2FoldChange <- results$logFC
results$minusLog10AdjustedPvalue <- -log10(results$FDR)

# Create a volcano plot with a significance threshold of FDR < 0.05
threshold <- 0.05
ggplot(results, aes(x = log2FoldChange, y = minusLog10AdjustedPvalue)) +
  geom_point(aes(color = FDR < threshold), alpha = 0.8, size = 2) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  theme(legend.position = "none") +
  xlab("log2(Fold Change)") +
  ylab("-log10(FDR)") +
  ggtitle("Volcano Plot")


