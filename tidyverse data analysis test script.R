library(airway)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(SummarizedExperiment)
library(EnhancedVolcano)
library(KEGGREST)


# -----------------------------
# 2. Load airway dataset
# -----------------------------
data("airway")
airway

# View sample metadata
colData(airway)

# -----------------------------
# 3. Create DESeq2 dataset
# -----------------------------
dds <- DESeqDataSet(airway, design = ~ dex)  # 'dex' = dexamethasone treatment

# Pre-filter: remove genes with very low counts
dds <- dds[rowSums(counts(dds)) > 1, ]

# -----------------------------
# 4. Run DESeq2
# -----------------------------
dds <- DESeq(dds)

# -----------------------------
# 5. Explore DESeq2 results
# -----------------------------
res <- results(dds)
summary(res)

# Order by p-value
resOrdered <- res[order(res$pvalue), ]
head(resOrdered)

# -----------------------------
# 6. Normalization & transformations
# -----------------------------
vsd <- vst(dds, blind = TRUE)    # Variance stabilizing transformation
rld <- rlog(dds, blind = TRUE)   # Regularized log transformation

# -----------------------------
# 7. Visualizations
# -----------------------------

## 7a. PCA plot
pca_data <- plotPCA(rld, intgroup = "dex", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
ggplot(pca_data, aes(PC1, PC2, color = dex)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of Airway Samples") +
  theme_minimal()

## 7b. MA plot
plotMA(res, main = "Airway DESeq2 MA Plot")

## 7c. Heatmap of top 20 variable genes
topGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
mat <- assay(vsd)[topGenes, ]
mat <- mat - rowMeans(mat)
pheatmap(mat, annotation_col = as.data.frame(colData(dds)[, "dex", drop = FALSE]),
         main = "Top 20 Variable Genes Heatmap",
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255))

# -----------------------------
# 8. Gene Set Enrichment Analysis (GSEA)
# -----------------------------
# Prepare named vector for GSEA (log2 fold change)
gene_list <- res$log2FoldChange
names(gene_list) <- rownames(res)
gene_list <- sort(na.omit(gene_list), decreasing = TRUE)

gse_res <- gseGO(geneList = gene_list,
                 OrgDb = org.Hs.eg.db,
                 ont = "BP",
                 keyType = "ENSEMBL",
                 eps = 1e-300)

# View top GSEA results
head(as.data.frame(gse_res))

# Dotplot of enriched pathways
dotplot(gse_res, showCategory = 10) + ggtitle("Top 10 Enriched GO Terms")

# -----------------------------
# 9. Export results
# -----------------------------
write.csv(as.data.frame(resOrdered), file = "Airway_DESeq2_results.csv")
write.csv(as.data.frame(gse_res), file = "Airway_GSEA_results.csv")

# -----------------------------
# 10. Volcano Plot of DE Genes
# -----------------------------
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                title = "Volcano plot of DE genes")

# -----------------------------
# 11. Sample Distance Heatmap
# -----------------------------
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         main = "Sample Distance Heatmap")

# -----------------------------
# 12. GO Overrepresentation Analysis (Upregulated Genes)
# -----------------------------

# Select significantly upregulated genes
sig_genes <- rownames(res)[which(res$padj < 0.05 & res$log2FoldChange > 0)]

# Enrichment in Biological Process (BP)
ego <- enrichGO(gene = sig_genes,
                OrgDb = org.Hs.eg.db,
                keyType = "ENSEMBL",
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05)

# Dotplot of top GO terms
dotplot(ego, showCategory = 10) + ggtitle("GO Enrichment - Upregulated Genes")


# -----------------------------
# 13. KEGG Pathway Enrichment
# -----------------------------

# Map ENSEMBL IDs to Entrez IDs
ensembl_genes_clean <- gsub("\\..*$", "", rownames(res))
entrez_map <- bitr(ensembl_genes_clean,
                   fromType = "ENSEMBL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

gene_list_entrez <- res$log2FoldChange[match(entrez_map$ENSEMBL, ensembl_genes_clean)]
names(gene_list_entrez) <- entrez_map$ENTREZID
gene_list_entrez <- sort(na.omit(gene_list_entrez), decreasing = TRUE)

# KEGG enrichment
kegg_res <- enrichKEGG(gene = names(gene_list_entrez),
                       organism = 'hsa',
                       keyType = 'ncbi-geneid',
                       pvalueCutoff = 0.05)

# Dotplot of top KEGG pathways
dotplot(kegg_res, showCategory = 10) + ggtitle("KEGG Pathway Enrichment")

# -----------------------------
# 14. Save individual plots as PNGs
# -----------------------------

# Volcano plot
ggsave("Volcano_plot.png", plot = last_plot(), width = 7, height = 5)

# Heatmap of top 20 variable genes (from step 7c)
pheatmap(mat,
         annotation_col = as.data.frame(colData(dds)[, "dex", drop = FALSE]),
         main = "Top 20 Variable Genes Heatmap",
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255),
         filename = "Top20_VariableGenes_Heatmap.png")

# GO enrichment dotplot
ggsave("GO_Enrichment_Dotplot.png", plot = last_plot(), width = 7, height = 5)

# KEGG enrichment dotplot
ggsave("KEGG_Enrichment_Dotplot.png", plot = last_plot(), width = 7, height = 5)

