library(airway)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(SummarizedExperiment)

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