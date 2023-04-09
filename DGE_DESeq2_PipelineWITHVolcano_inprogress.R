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
# -----------------------------------------------------------

#MAKING DATAFRAMES FOR DOWNSTREAM USE
#-------
# read in the data ---------
dat <- read.table(file = "data/bulk/GSE174465_Normalized_counts.txt", 
                  sep = "\t", header = T, stringsAsFactors = FALSE)
# get metadata --------
gse <- getGEO(GEO = 'GSE174465', GSEMatrix = T)

metadata <- pData(phenoData(gse[[1]]))
head(metadata)

# select, mutate, rename ------------
metadata.modified <- metadata %>%
  select(1,8,44,45,46) %>%
  rename_with(.cols = 2, ~"cell_type") %>%
  rename_with(.cols = 3, ~"clone") %>%
  rename_with(.cols = 4, ~"cell_info") %>%
  rename_with(.cols = 5, ~"sex") %>%
  mutate(cell_type = gsub("human ", "", cell_type))
head(metadata.modified)

# looking at gene expression data ---------
head(dat)

# reshaping data - from wide to long--------
dat.long <- dat %>%
  rename_with(.cols = 1, ~"gene") %>%
  gather(key = 'samples', value = 'norm.counts', -gene)
head(dat.long)

# join dataframes = dat.long + metadata.modified
dat.long <- dat.long %>%
  left_join(., metadata.modified, by = c("samples" = "title")) 

#------------
#DGE updates
#-----------
# Step 1: preparing count data ----------------

#updating counts data
head(dat)
counts_data <- dat
rownames(counts_data) <- counts_data$X
counts_data$X <- NULL
head(counts_data)

#updating design matrix
colData <- metadata.modified
rownames(colData) <- colData$title
colData$title <- NULL
head(colData)

# making sure the row names in colData matches to column names in counts_data
all(colnames(counts_data) %in% rownames(colData))

# are they in the same order?
all(colnames(counts_data) == rownames(colData))


# Step 2: construct a DESeqDataSet object ----------

dds <- DESeqDataSetFromMatrix(countData = round(counts_data),
                              colData = colData,
                              design = ~ cell_type)
dds <- DESeq(dds)

res <- results(dds, contrast = c("cell_type", "PGCLC", "EG"))

res <- data.frame(log2FoldChange = res$log2FoldChange,
                  lfcSE = res$lfcSE,
                  pvalue = res$pvalue,
                  padj = res$padj,
                  row.names = rownames(res))

#Make a set of DEGs to create a .csv file 
# Subset the differentially expressed genes
DE_genes <- res[abs(res$log2FoldChange) > 1 & res$padj < 0.05, ]

# Add a column for the sign of the log2 fold change
DE_genes$log2FC_sign <- ifelse(DE_genes$log2FoldChange > 0, "up", "down")

# Sort the genes by log2 fold change absolute value
DE_genes <- DE_genes[order(-abs(DE_genes$log2FoldChange)), ]

# Export the differentially expressed genes to a CSV file
write.csv(DE_genes, "bulk_analysis/results/DE_genes_info.csv", row.names = TRUE)


EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                legendPosition = 'right',
                title = 'Differentially expressed genes: EGCLC vs PGCLC')

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = c('SOX17','SOX2','NANOS3',
                              'PRDM1', 'TFAP2C','PIWIL2', 'DPPA3',
                              'DLK1','CD24'),
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'Key Genes of Interest',
                FCcutoff = 2.0,
                pointSize = 4.0,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')


EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = c('SSEA-4'),
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'Key Genes of Interest')

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'Key Genes of Interest',
                FCcutoff = 2.0,
                pointSize = 4.0,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')

#-------------- in prog
# extract only up-regulated genes
up_genes <- subset(DE_genes, DE_genes$log2FC_sign == "up")

# convert gene symbols to Entrez IDs using org.Hs.eg.db database
up_genes_entrez <- bitr(rownames(up_genes), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

# perform GO analysis on up-regulated genes
go_enrichment_upgenes <- enrichGO(up_genes_entrez$ENTREZID, 
                                  OrgDb=org.Hs.eg.db, 
                                  keyType="ENTREZID", 
                                  ont="BP", 
                                  pvalueCutoff=0.05, 
                                  qvalueCutoff=0.1)

head(go_enrichment_upgenes)

barplot(go_enrichment_upgenes,
        showCategory=20, # show the top 10 enriched GO terms
        title="Enriched Biological Processes (PGCLCs)")


# extract only up-regulated genes
down_genes <- subset(DE_genes, DE_genes$log2FC_sign == "down")

# convert gene symbols to Entrez IDs using org.Hs.eg.db database
down_genes_entrez <- bitr(rownames(down_genes), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

# perform GO analysis on up-regulated genes
go_enrichment_downgenes <- enrichGO(down_genes_entrez$ENTREZID, 
                                    OrgDb=org.Hs.eg.db, 
                                    keyType="ENTREZID", 
                                    ont="BP", 
                                    pvalueCutoff=0.05, 
                                    qvalueCutoff=0.1)
head(go_enrichment_downgenes)
barplot(go_enrichment_downgenes,
        showCategory=20, # show the top 10 enriched GO terms
        title="Enriched Biological Processes (EGCLCs)")


# perform GO analysis on up-regulated genes
go_enrichment_downgenes <- enrichGO(down_genes_entrez$ENTREZID, 
                                    OrgDb=org.Hs.eg.db, 
                                    keyType="ENTREZID", 
                                    ont="BP", 
                                    pvalueCutoff=0.05, 
                                    qvalueCutoff=0.1)




 edox2 <- pairwise_termsim(go_enrichment_downgenes)
 p1 <- treeplot(edox2)
 p2 <- treeplot(edox2, hclust_method = "average")
 aplot::plot_list(p1, p2, tag_levels='A')
 
 
 
 #----------
 
 go_enrichment_downgenes_all <- enrichGO(down_genes_entrez$ENTREZID, 
                                     OrgDb=org.Hs.eg.db, 
                                     keyType="ENTREZID", 
                                     ont="ALL", 
                                     pvalueCutoff=0.05, 
                                     qvalueCutoff=0.1)
 barplot(go_enrichment_downgenes_all,
         showCategory=20, # show the top 10 enriched GO terms
         title="Enriched Processes (EGCLCs)")
 edox2 <- pairwise_termsim(go_enrichment_upgenes_all)
 p1 <- treeplot(edox2)
 p2 <- treeplot(edox2, hclust_method = "average")
 aplot::plot_list(p1, p2, tag_levels='A')
 
 
 go_enrichment_upgenes_all <- enrichGO(up_genes_entrez$ENTREZID, 
                                         OrgDb=org.Hs.eg.db, 
                                         keyType="ENTREZID", 
                                         ont="ALL", 
                                         pvalueCutoff=0.05, 
                                         qvalueCutoff=0.1)
 
 
 
 