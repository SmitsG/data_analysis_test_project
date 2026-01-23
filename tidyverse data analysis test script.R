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
