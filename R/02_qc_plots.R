# R/02_qc_plots.R
suppressPackageStartupMessages({
  library(ggplot2)
  library(DESeq2)
  library(SummarizedExperiment)
  library(pheatmap)
  library(RColorBrewer)
  library(vsn)
  library(matrixStats)
})

plot_pca <- function(rld, intgroup = "dex", title = "PCA (rlog)") {
  pca_data <- DESeq2::plotPCA(rld, intgroup = intgroup, returnData = TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"))
  
  ggplot(pca_data, aes(PC1, PC2, color = .data[[intgroup]])) +
    geom_point(size = 4) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle(title) +
    theme_minimal()
}

plot_ma_png <- function(res, path) {
  save_base_png(path, {
    DESeq2::plotMA(res, main = "MA plot")
  }, width = 1100, height = 850)
}

plot_dispersion_png <- function(dds, path) {
  save_base_png(path, {
    DESeq2::plotDispEsts(dds, main = "Dispersion estimates")
  }, width = 1100, height = 850)
}

plot_cooks_png <- function(dds, path) {
  save_base_png(path, {
    boxplot(log10(SummarizedExperiment::assays(dds)[["cooks"]]),
            main = "Log10 Cook's distances", las = 2)
  }, width = 1200, height = 850)
}

plot_mean_sd <- function(vsd, rld) {
  list(
    vsd = vsn::meanSdPlot(SummarizedExperiment::assay(vsd))$gg + ggtitle("Mean-SD (VST)"),
    rld = vsn::meanSdPlot(SummarizedExperiment::assay(rld))$gg + ggtitle("Mean-SD (rlog)")
  )
}

plot_sample_distance_heatmap <- function(vsd, annotation_col, path) {
  sampleDists <- stats::dist(t(SummarizedExperiment::assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  cols <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)
  
  pheatmap::pheatmap(
    sampleDistMatrix,
    col = cols,
    annotation_col = annotation_col,
    main = "Sample-to-sample distance (VST)",
    filename = path
  )
  path
}

plot_top_variable_heatmap <- function(vsd, annotation_col, n = 20, path) {
  topGenes <- head(order(matrixStats::rowVars(SummarizedExperiment::assay(vsd)), decreasing = TRUE), n)
  mat <- SummarizedExperiment::assay(vsd)[topGenes, , drop = FALSE]
  mat <- mat - rowMeans(mat)
  
  pheatmap::pheatmap(
    mat,
    annotation_col = annotation_col,
    main = paste0("Top ", n, " variable genes (VST)"),
    color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(255),
    filename = path
  )
  path
}
