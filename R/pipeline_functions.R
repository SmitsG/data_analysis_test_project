# R/pipeline_functions.R
suppressPackageStartupMessages({
  library(here)
})

source(here::here("R", "00_utils.R"))
source(here::here("R", "01_deseq2.R"))
source(here::here("R", "02_qc_plots.R"))
source(here::here("R", "03_enrichment.R"))
source(here::here("R", "04_reporting.R"))

suppressPackageStartupMessages({
  library(airway)
  library(DESeq2)
  library(SummarizedExperiment)
  library(ggplot2)
  library(EnhancedVolcano)
})

run_airway_pipeline <- function(alpha = 0.05,
                                padj_cutoff = 0.05,
                                shrink = TRUE,
                                coef = "dex_trt_vs_untrt") {
  dirs <- init_output_dirs()
  log_path <- file.path(dirs$results, "pipeline.log")
  safe_remove(log_path)
  
  log_line(log_path, "Starting airway pipeline")
  
  data("airway", package = "airway")
  se <- airway
  log_line(log_path, "Loaded airway dataset | samples:", ncol(se), "genes:", nrow(se))
  
  dds <- DESeq2::DESeqDataSet(se, design = ~ dex)
  dds <- prefilter_dds(dds, min_total = 1)
  dds <- run_deseq2(dds)
  log_line(log_path, "DESeq2 fit complete | genes:", nrow(dds), "samples:", ncol(dds))
  
  res <- get_deseq_results(dds, alpha = alpha)
  save_csv(as.data.frame(res), file.path(dirs$results, "deseq2_results.csv"))
  
  res_shrunk <- if (isTRUE(shrink)) shrink_lfc_safe(dds, coef = coef) else res
  res_shrunk <- res_shrunk[order(res_shrunk$pvalue), ]
  save_csv(as.data.frame(res_shrunk), file.path(dirs$results, "deseq2_results_shrunk.csv"))
  
  tx <- make_transforms(dds)
  vsd <- tx$vsd
  rld <- tx$rld
  ann <- as.data.frame(SummarizedExperiment::colData(dds)[, "dex", drop = FALSE])
  
  pca <- plot_pca(rld, intgroup = "dex", title = "PCA (rlog) - Airway")
  save_gg(pca, file.path(dirs$plots, "PCA_plot.png"))
  
  plot_ma_png(res, file.path(dirs$plots, "MA_plot.png"))
  plot_dispersion_png(dds, file.path(dirs$plots, "Dispersion_plot.png"))
  plot_cooks_png(dds, file.path(dirs$plots, "Cooks_distance_boxplot.png"))
  plot_sample_distance_heatmap(vsd, ann, file.path(dirs$plots, "Sample_distance_heatmap.png"))
  plot_top_variable_heatmap(vsd, ann, n = 20, path = file.path(dirs$plots, "Top20_VariableGenes_Heatmap.png"))
  
  ms <- plot_mean_sd(vsd, rld)
  save_gg(ms$vsd, file.path(dirs$plots, "MeanSD_VST.png"))
  save_gg(ms$rld, file.path(dirs$plots, "MeanSD_RLD.png"))
  
  volcano <- EnhancedVolcano::EnhancedVolcano(
    res_shrunk,
    lab = rownames(res_shrunk),
    x = "log2FoldChange",
    y = "pvalue",
    pCutoff = 0.05,
    FCcutoff = 1.5,
    title = "Volcano plot - Airway (shrunk LFC if available)"
  )
  save_gg(volcano, file.path(dirs$plots, "Volcano_plot.png"), width = 7, height = 5)
  
  sig <- get_sig_genes(res, padj_cutoff = padj_cutoff)
  save_csv(data.frame(ENSEMBL = sig$up), file.path(dirs$results, "Upregulated_genes.csv"))
  save_csv(data.frame(ENSEMBL = sig$down), file.path(dirs$results, "Downregulated_genes.csv"))
  log_line(log_path, "Sig genes | up:", length(sig$up), "down:", length(sig$down))
  
  ego <- run_go_ora(sig$up, ont = "BP")
  save_csv(as.data.frame(ego), file.path(dirs$results, "GO_ORA_upregulated.csv"))
  
  go_dot <- safe_dotplot(ego, "GO BP ORA - Upregulated", showCategory = 10)
  if (!is.null(go_dot)) save_gg(go_dot, file.path(dirs$plots, "GO_Enrichment_Dotplot.png"), width = 8, height = 6)
  
  go_emap <- go_emapplot_safe(ego, showCategory = 10)
  if (!is.null(go_emap)) save_gg(go_emap, file.path(dirs$plots, "GO_emapplot.png"), width = 10, height = 8)
  
  go_cnet <- go_cnetplot_safe(ego, showCategory = 5)
  if (!is.null(go_cnet)) save_gg(go_cnet, file.path(dirs$plots, "GO_cnetplot.png"), width = 10, height = 8)
  
  rank_obj <- build_entrez_rank(res_shrunk)
  lfc_rank <- rank_obj$lfc_rank
  symbols <- add_symbols(names(lfc_rank))
  
  ranked_df <- data.frame(
    EntrezID = names(lfc_rank),
    Symbol = unname(symbols),
    log2FC = as.numeric(lfc_rank)
  )
  save_csv(ranked_df, file.path(dirs$results, "Ranked_Genes_for_GSEA_with_SYMBOL.csv"))
  log_line(log_path, "Ranked Entrez genes:", nrow(ranked_df))
  
  kegg_ora <- run_kegg_ora(names(lfc_rank), organism = "hsa")
  save_csv(as.data.frame(kegg_ora), file.path(dirs$results, "KEGG_ORA.csv"))
  
  kegg_dot <- safe_dotplot(kegg_ora, "KEGG ORA", showCategory = 10)
  if (!is.null(kegg_dot)) save_gg(kegg_dot, file.path(dirs$plots, "KEGG_Enrichment_Dotplot.png"), width = 8, height = 6)
  
  gse_go <- run_go_gsea(lfc_rank, ont = "BP")
  save_csv(as.data.frame(gse_go), file.path(dirs$results, "GO_GSEA.csv"))
  
  gse_kegg <- run_kegg_gsea(lfc_rank, organism = "hsa")
  save_csv(as.data.frame(gse_kegg), file.path(dirs$results, "KEGG_GSEA.csv"))
  
  ridge_go <- safe_ridgeplot(gse_go, "GSEA GO ridgeplot", showCategory = 15)
  ridge_kegg <- safe_ridgeplot(gse_kegg, "GSEA KEGG ridgeplot", showCategory = 15)
  if (!is.null(ridge_go)) save_gg(ridge_go, file.path(dirs$plots, "GSEA_GO_ridgeplot.png"), width = 9, height = 6)
  if (!is.null(ridge_kegg)) save_gg(ridge_kegg, file.path(dirs$plots, "GSEA_KEGG_ridgeplot.png"), width = 9, height = 6)
  
  dot_go_sign <- safe_dotplot_split(gse_go, "GSEA GO dotplot (split sign)", showCategory = 10)
  dot_kegg_sign <- safe_dotplot_split(gse_kegg, "GSEA KEGG dotplot (split sign)", showCategory = 10)
  if (!is.null(dot_go_sign)) save_gg(dot_go_sign, file.path(dirs$plots, "GSEA_GO_dotplot_sign.png"), width = 10, height = 6)
  if (!is.null(dot_kegg_sign)) save_gg(dot_kegg_sign, file.path(dirs$plots, "GSEA_KEGG_dotplot_sign.png"), width = 10, height = 6)
  
  pdf_path <- file.path(dirs$plots, "Airway_DE_Analysis_All_Plots.pdf")
  write_pdf_report(
    pdf_path,
    ggplots = list(
      pca, volcano, go_dot, kegg_dot,
      ms$vsd, ms$rld,
      ridge_go, ridge_kegg,
      dot_go_sign, dot_kegg_sign,
      go_emap, go_cnet
    ),
    base_plot_fns = list(function() DESeq2::plotMA(res, main = "MA plot - Airway"))
  )
  
  log_line(log_path, "Saved PDF report:", pdf_path)
  log_line(log_path, "Finished airway pipeline")
  
  list(
    dds = dds,
    res = res,
    res_shrunk = res_shrunk,
    vsd = vsd,
    rld = rld,
    ego = ego,
    kegg_ora = kegg_ora,
    gse_go = gse_go,
    gse_kegg = gse_kegg,
    ranked = ranked_df,
    dirs = dirs,
    log = log_path
  )
}
