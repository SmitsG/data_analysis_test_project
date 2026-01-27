# ==============================================
# Airway DESeq2 + Enrichment Pipeline (Improved)
# ==============================================

# -----------------------------
# 0) Setup: packages + folders
# -----------------------------
suppressPackageStartupMessages({
  library(airway)
  library(DESeq2)
  library(SummarizedExperiment)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(enrichplot)
  library(DOSE)
  library(EnhancedVolcano)
  library(vsn)
  library(matrixStats)
  library(dplyr)
  library(here)
})

# Create output dirs (safe if already exist)
dir.create(here("plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("results"), showWarnings = FALSE, recursive = TRUE)

# Helper: safe ggsave
save_plot <- function(plot_obj, filename, w = 7, h = 5) {
  ggplot2::ggsave(here("plots", filename), plot = plot_obj, width = w, height = h)
  invisible(here("plots", filename))
}

# Helper: save base plot to png
save_png <- function(filename, expr, width = 1000, height = 800) {
  grDevices::png(here("plots", filename), width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
  invisible(here("plots", filename))
}

# Helper: assert a condition with readable message
assert_that <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

# -----------------------------
# 1) Load dataset + inspect
# -----------------------------
data("airway")
message("Samples: ", ncol(airway), " | Genes: ", nrow(airway))
print(colData(airway))

# -----------------------------
# 2) Build DESeq2 dataset
# -----------------------------
dds <- DESeqDataSet(airway, design = ~ dex)

# Prefilter: remove very low count genes
dds <- dds[rowSums(counts(dds)) > 1, ]

# Run DESeq2
dds <- DESeq(dds)

# Results (default: dex treated vs untreated)
res <- results(dds)
resOrdered <- res[order(res$pvalue), ]

# Save DE results table
write.csv(as.data.frame(resOrdered), here("results", "Airway_DESeq2_results.csv"))

# Optional: LFC shrinkage (recommended for nicer ranking/plots)
# Uses DESeq2's built-in shrinkers; type="apeglm" requires apeglm installed.
res_shrunk <- tryCatch(
  lfcShrink(dds, coef = "dex_trt_vs_untrt", type = "apeglm"),
  error = function(e) {
    message("lfcShrink(apeglm) not available (install apeglm). Using normal results.")
    res
  }
)
write.csv(as.data.frame(res_shrunk[order(res_shrunk$pvalue), ]), here("results", "Airway_DESeq2_results_shrunk.csv"))

# -----------------------------
# 3) Transformations for QC/plots
# -----------------------------
vsd <- vst(dds, blind = TRUE)
rld <- rlog(dds, blind = TRUE)

# -----------------------------
# 4) QC plots
# -----------------------------

# 4a) PCA
pca_data <- plotPCA(rld, intgroup = "dex", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = dex)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA (rlog) - Airway") +
  theme_minimal()
save_plot(pca_plot, "PCA_plot.png")

# 4b) MA plot (base plot -> png)
save_png("MA_plot.png", plotMA(res, main = "MA plot - Airway"))

# 4c) Sample-to-sample distance heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
dist_cols <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(
  sampleDistMatrix,
  col = dist_cols,
  main = "Sample-to-sample distance (VST)",
  filename = here("plots", "Sample_distance_heatmap.png")
)

# 4d) Mean-SD plots
ms_vsd <- meanSdPlot(assay(vsd))$gg + ggtitle("Mean-SD (VST)")
ms_rld <- meanSdPlot(assay(rld))$gg + ggtitle("Mean-SD (rlog)")
save_plot(ms_vsd, "MeanSD_VST.png")
save_plot(ms_rld, "MeanSD_RLD.png")

# 4e) Cook's distance (influential counts)
save_png("Cooks_distance_boxplot.png", {
  boxplot(log10(assays(dds)[["cooks"]]), main = "Log10 Cook's distances", las = 2)
})

# 4f) Dispersion plot
save_png("Dispersion_plot.png", plotDispEsts(dds, main = "Dispersion estimates"))

# 4g) Top variable genes heatmap (top 20)
topGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
mat_top <- assay(vsd)[topGenes, , drop = FALSE]
mat_top <- mat_top - rowMeans(mat_top)
ann <- as.data.frame(colData(dds)[, "dex", drop = FALSE])
pheatmap(
  mat_top,
  annotation_col = ann,
  main = "Top 20 variable genes (VST)",
  color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255),
  filename = here("plots", "Top20_VariableGenes_Heatmap.png")
)

# -----------------------------
# 5) Volcano plot (use shrunk LFC if available)
# -----------------------------
volcano <- EnhancedVolcano(
  res_shrunk,
  lab = rownames(res_shrunk),
  x = "log2FoldChange",
  y = "pvalue",
  pCutoff = 0.05,
  FCcutoff = 1.5,
  title = "Volcano plot - Airway (LFC shrink if available)"
)
save_plot(volcano, "Volcano_plot.png", w = 7, h = 5)

# -----------------------------
# 6) Define significant gene sets (up/down)
# -----------------------------
sig_up <- rownames(res)[which(res$padj < 0.05 & res$log2FoldChange > 0)]
sig_down <- rownames(res)[which(res$padj < 0.05 & res$log2FoldChange < 0)]

write.csv(data.frame(ENSEMBL = sig_up), here("results", "Upregulated_genes.csv"), row.names = FALSE)
write.csv(data.frame(ENSEMBL = sig_down), here("results", "Downregulated_genes.csv"), row.names = FALSE)

# -----------------------------
# 7) GO ORA (upregulated)
# -----------------------------
ego <- enrichGO(
  gene = sig_up,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
write.csv(as.data.frame(ego), here("results", "GO_ORA_upregulated.csv"), row.names = FALSE)

go_dot <- dotplot(ego, showCategory = 10) + ggtitle("GO BP ORA - Upregulated")
save_plot(go_dot, "GO_Enrichment_Dotplot.png", w = 8, h = 6)

# GO term network + emap
ego_simplified <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
ego_sim <- pairwise_termsim(ego_simplified)

go_emap <- emapplot(ego_sim, showCategory = 10) + ggtitle("GO functional clusters (emapplot)")
save_plot(go_emap, "GO_emapplot.png", w = 10, h = 8)

go_cnet <- cnetplot(ego, showCategory = 5) + ggtitle("GO cnetplot (top 5)")
save_plot(go_cnet, "GO_cnetplot.png", w = 10, h = 8)

# -----------------------------
# 8) Map ENSEMBL -> ENTREZ and run KEGG ORA + GSEA
# -----------------------------
ensembl_clean <- gsub("\\..*$", "", rownames(res))
entrez_map <- bitr(
  ensembl_clean,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Build ranked Entrez list using shrunk LFC if possible
lfc_rank <- res_shrunk$log2FoldChange[match(entrez_map$ENSEMBL, ensembl_clean)]
names(lfc_rank) <- entrez_map$ENTREZID

# Clean ranked list: drop NA, drop duplicates, sort decreasing
lfc_rank <- lfc_rank[!is.na(lfc_rank)]
lfc_rank <- lfc_rank[!duplicated(names(lfc_rank))]
lfc_rank <- sort(lfc_rank, decreasing = TRUE)

assert_that(length(lfc_rank) > 1000, "Ranked Entrez list is unexpectedly small; mapping failed?")

# Save ranked list (with SYMBOL column)
symbols <- mapIds(org.Hs.eg.db, keys = names(lfc_rank), column = "SYMBOL",
                  keytype = "ENTREZID", multiVals = "first")
ranked_df <- data.frame(EntrezID = names(lfc_rank), Symbol = unname(symbols), log2FC = as.numeric(lfc_rank))
write.csv(ranked_df, here("results", "Ranked_Genes_for_GSEA_with_SYMBOL.csv"), row.names = FALSE)

# KEGG ORA (uses Entrez IDs)
kegg_ora <- enrichKEGG(
  gene = names(lfc_rank),
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05
)
write.csv(as.data.frame(kegg_ora), here("results", "KEGG_ORA.csv"), row.names = FALSE)

kegg_dot <- dotplot(kegg_ora, showCategory = 10) + ggtitle("KEGG ORA")
save_plot(kegg_dot, "KEGG_Enrichment_Dotplot.png", w = 8, h = 6)

# GSEA GO (Entrez list)
gse_go <- gseGO(
  geneList = lfc_rank,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "ENTREZID",
  eps = 1e-300
)
write.csv(as.data.frame(gse_go), here("results", "GO_GSEA.csv"), row.names = FALSE)

# GSEA KEGG
gse_kegg <- gseKEGG(
  geneList = lfc_rank,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500
)
write.csv(as.data.frame(gse_kegg), here("results", "KEGG_GSEA.csv"), row.names = FALSE)

# Ridgeplots + dotplots (split sign)
ridge_go <- ridgeplot(gse_go, showCategory = 15) + ggtitle("GSEA GO ridgeplot")
ridge_kegg <- ridgeplot(gse_kegg, showCategory = 15) + ggtitle("GSEA KEGG ridgeplot")
save_plot(ridge_go, "GSEA_GO_ridgeplot.png", w = 9, h = 6)
save_plot(ridge_kegg, "GSEA_KEGG_ridgeplot.png", w = 9, h = 6)

dot_go_sign <- dotplot(gse_go, showCategory = 10, split = ".sign") + DOSE::facet_grid(.~.sign)
dot_kegg_sign <- dotplot(gse_kegg, showCategory = 10, split = ".sign") + DOSE::facet_grid(.~.sign)
save_plot(dot_go_sign, "GSEA_GO_dotplot_sign.png", w = 10, h = 6)
save_plot(dot_kegg_sign, "GSEA_KEGG_dotplot_sign.png", w = 10, h = 6)

# -----------------------------
# 9) One PDF report with all key plots
# -----------------------------
pdf(here("plots", "Airway_DE_Analysis_All_Plots.pdf"), width = 10, height = 8)

print(pca_plot)
plotMA(res, main = "MA plot - Airway")

print(volcano)

pheatmap(mat_top,
         annotation_col = ann,
         main = "Top 20 variable genes (VST)",
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255))

pheatmap(sampleDistMatrix, col = dist_cols, main = "Sample-to-sample distance (VST)")

print(go_dot)
print(kegg_dot)
print(ms_vsd); print(ms_rld)
print(ridge_go); print(ridge_kegg)
print(dot_go_sign); print(dot_kegg_sign)
print(go_emap)
print(go_cnet)

dev.off()

message("Done. Outputs saved to: ", here("plots"), " and ", here("results"))
