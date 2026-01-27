# R/functions.R

# ---- Packages used inside functions ----
suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(DOSE)
  library(EnhancedVolcano)
  library(vsn)
  library(dplyr)
})

# ---- Utilities ----
ensure_dirs <- function() {
  dirs <- c(here::here("plots"), here::here("results"))
  for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  invisible(dirs)
}

save_csv <- function(df, path) {
  utils::write.csv(df, file = path, row.names = FALSE)
  path
}

# Save ggplot-like object safely
save_gg <- function(p, path, width = 7, height = 5) {
  ggplot2::ggsave(filename = path, plot = p, width = width, height = height)
  path
}

# Save a base R plot call to PNG
save_base_png <- function(path, expr, width = 800, height = 600) {
  grDevices::png(path, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
  path
}

# ---- DESeq2 core ----
make_dds <- function(se, design = ~ dex, min_count_sum = 1) {
  dds <- DESeqDataSet(se, design = design)
  dds <- dds[rowSums(counts(dds)) > min_count_sum, ]
  dds
}

run_deseq <- function(dds) {
  DESeq(dds)
}

get_results <- function(dds) {
  res <- results(dds)
  resOrdered <- res[order(res$pvalue), ]
  list(res = res, resOrdered = resOrdered)
}

make_transforms <- function(dds) {
  list(
    vsd = vst(dds, blind = TRUE),
    rld = rlog(dds, blind = TRUE)
  )
}

# ---- Plots ----
plot_pca <- function(rld, intgroup = "dex") {
  pca_data <- plotPCA(rld, intgroup = intgroup, returnData = TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"))
  
  ggplot(pca_data, aes(PC1, PC2, color = .data[[intgroup]])) +
    geom_point(size = 4) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle("PCA of Samples") +
    theme_minimal()
}

plot_ma_png <- function(res, path) {
  save_base_png(path, {
    plotMA(res, main = "MA Plot")
  }, width = 900, height = 700)
}

heatmap_top_var <- function(vsd, dds, n = 20, path) {
  topGenes <- head(order(matrixStats::rowVars(assay(vsd)), decreasing = TRUE), n)
  mat <- assay(vsd)[topGenes, , drop = FALSE]
  mat <- mat - rowMeans(mat)
  
  ann <- as.data.frame(colData(dds)[, "dex", drop = FALSE])
  
  pheatmap::pheatmap(
    mat,
    annotation_col = ann,
    main = paste0("Top ", n, " Variable Genes Heatmap"),
    color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255),
    filename = path
  )
  path
}

plot_volcano <- function(res) {
  EnhancedVolcano::EnhancedVolcano(
    res,
    lab = rownames(res),
    x = "log2FoldChange",
    y = "pvalue",
    pCutoff = 0.05,
    FCcutoff = 1.5,
    title = "Volcano plot of DE genes"
  )
}

sample_distance_heatmap <- function(vsd, path) {
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  
  pheatmap::pheatmap(
    sampleDistMatrix,
    col = colors,
    main = "Sample-to-sample distance heatmap",
    filename = path
  )
  path
}

# ---- ORA (GO) ----
get_sig_genes <- function(res, padj_cut = 0.05) {
  list(
    up = rownames(res)[which(res$padj < padj_cut & res$log2FoldChange > 0)],
    down = rownames(res)[which(res$padj < padj_cut & res$log2FoldChange < 0)]
  )
}

ora_go_up <- function(sig_up) {
  enrichGO(
    gene = sig_up,
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  )
}

# ---- ID mapping + KEGG ----
map_ensembl_to_entrez <- function(res) {
  ensembl_clean <- gsub("\\..*$", "", rownames(res))
  m <- bitr(
    ensembl_clean,
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  # match log2FC to mapped ENSEMBL
  lfc <- res$log2FoldChange[match(m$ENSEMBL, ensembl_clean)]
  names(lfc) <- m$ENTREZID
  
  # remove NAs & duplicates, ensure decreasing sort
  lfc <- lfc[!is.na(lfc)]
  lfc <- lfc[!duplicated(names(lfc))]
  lfc <- sort(lfc, decreasing = TRUE)
  
  list(entrez_map = m, gene_list_entrez = lfc)
}

ora_kegg <- function(gene_list_entrez) {
  enrichKEGG(
    gene = names(gene_list_entrez),
    organism = "hsa",
    keyType = "ncbi-geneid",
    pvalueCutoff = 0.05
  )
}

# ---- GSEA ----
gsea_go <- function(gene_list_entrez) {
  gseGO(
    geneList = gene_list_entrez,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    keyType = "ENTREZID",
    eps = 1e-300
  )
}

gsea_kegg <- function(gene_list_entrez) {
  gseKEGG(
    geneList = gene_list_entrez,
    organism = "hsa",
    keyType = "ncbi-geneid",
    pvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500
  )
}

# ---- Extra QC ----
mean_sd_plots <- function(vsd, rld, dds) {
  list(
    vsd = meanSdPlot(assay(vsd))$gg + ggtitle("Mean-SD Plot (VST)"),
    rld = meanSdPlot(assay(rld))$gg + ggtitle("Mean-SD Plot (RLD)"),
    ntd = meanSdPlot(assay(dds))$gg + ggtitle("Mean-SD Plot (untransformed)")
  )
}

cooks_png <- function(dds, path) {
  save_base_png(path, {
    boxplot(log10(assays(dds)[["cooks"]]), main = "Log10 Cook's Distance", las = 2)
  }, width = 1100, height = 700)
}

dispersion_png <- function(dds, path) {
  save_base_png(path, {
    plotDispEsts(dds, main = "Dispersion Estimates")
  }, width = 1100, height = 700)
}

ma_multipanel_png <- function(dds, path) {
  resGA <- results(dds, lfcThreshold = 0.5, altHypothesis = "greaterAbs")
  resLA <- results(dds, lfcThreshold = 0.5, altHypothesis = "lessAbs")
  resG  <- results(dds, lfcThreshold = 0.5, altHypothesis = "greater")
  resL  <- results(dds, lfcThreshold = 0.5, altHypothesis = "less")
  
  save_base_png(path, {
    par(mfrow = c(2,2), mar = c(2,2,1,1))
    drawLines <- function() abline(h = c(-0.5, 0.5), col = "dodgerblue", lwd = 2)
    plotMA(resGA); drawLines()
    plotMA(resLA); drawLines()
    plotMA(resG);  drawLines()
    plotMA(resL);  drawLines()
  }, width = 1400, height = 900)
}

# ---- Enrichment maps / networks ----
go_cnetplot <- function(ego, gene_list_entrez) {
  # color by fold change where possible (use named vector)
  enrichplot::cnetplot(
    ego,
    showCategory = 5,
    color.params = list(foldChange = gene_list_entrez)
  ) + ggtitle("GO Enrichment Network")
}

go_emapplot <- function(ego) {
  ego_simplified <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
  ego_sim <- pairwise_termsim(ego_simplified)
  enrichplot::emapplot(ego_sim, showCategory = 10) + ggtitle("GO Term Functional Clusters")
}

ora_full_emap <- function(res) {
  genes <- rownames(res)[which(res$padj < 0.05)]
  ego_full <- enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = "BP"
  )
  ego_sim <- pairwise_termsim(ego_full)
  enrichplot::emapplot(ego_sim, showCategory = 10) + ggtitle("Full ORA Enrichment Map")
}

# ---- Save ranked gene list with SYMBOL ----
save_ranked_genes <- function(gene_list_entrez, path_basic, path_symbols) {
  df_basic <- data.frame(EntrezID = names(gene_list_entrez), log2FC = as.numeric(gene_list_entrez))
  save_csv(df_basic, path_basic)
  
  symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = names(gene_list_entrez),
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  )
  df_sym <- data.frame(
    EntrezID = names(gene_list_entrez),
    Symbol = unname(symbols),
    log2FC = as.numeric(gene_list_entrez)
  )
  save_csv(df_sym, path_symbols)
  
  list(basic = path_basic, symbols = path_symbols)
}

# ---- PDF report ----
save_pdf_report <- function(path, res, dds, vsd, ego, kegg_res,
                            mean_sd, ridge_go, ridge_kegg, dot_go, dot_kegg, emap_full) {
  grDevices::pdf(path, width = 10, height = 8)
  on.exit(grDevices::dev.off(), add = TRUE)
  
  # EnhancedVolcano returns a ggplot-like object
  print(plot_volcano(res))
  
  # pheatmap draws directly
  topGenes <- head(order(matrixStats::rowVars(assay(vsd)), decreasing = TRUE), 20)
  mat <- assay(vsd)[topGenes, , drop = FALSE]
  mat <- mat - rowMeans(mat)
  ann <- as.data.frame(colData(dds)[, "dex", drop = FALSE])
  pheatmap::pheatmap(
    mat,
    annotation_col = ann,
    main = "Top 20 Variable Genes Heatmap",
    color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)
  )
  
  print(dotplot(ego, showCategory = 10) + ggtitle("GO Enrichment - Upregulated Genes"))
  print(dotplot(kegg_res, showCategory = 10) + ggtitle("KEGG Pathway Enrichment"))
  
  print(mean_sd$vsd); print(mean_sd$rld); print(mean_sd$ntd)
  print(ridge_go); print(ridge_kegg)
  print(dot_go); print(dot_kegg)
  print(emap_full)
  
  path
}
