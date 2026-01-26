# R/03_enrichment.R
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(AnnotationDbi)
  library(enrichplot)
  library(DOSE)
  library(ggplot2)
})

get_sig_genes <- function(res, padj_cutoff = 0.05) {
  up <- rownames(res)[which(res$padj < padj_cutoff & res$log2FoldChange > 0)]
  down <- rownames(res)[which(res$padj < padj_cutoff & res$log2FoldChange < 0)]
  list(up = up, down = down)
}

run_go_ora <- function(genes_ensembl,
                       ont = "BP",
                       p_adj_method = "BH",
                       qvalue_cutoff = 0.05) {
  clusterProfiler::enrichGO(
    gene = genes_ensembl,
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = ont,
    pAdjustMethod = p_adj_method,
    qvalueCutoff = qvalue_cutoff
  )
}

build_entrez_rank <- function(res_for_rank) {
  ensembl_clean <- gsub("\\..*$", "", rownames(res_for_rank))
  
  entrez_map <- clusterProfiler::bitr(
    ensembl_clean,
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  )
  
  lfc <- res_for_rank$log2FoldChange[match(entrez_map$ENSEMBL, ensembl_clean)]
  names(lfc) <- entrez_map$ENTREZID
  
  lfc <- lfc[!is.na(lfc)]
  lfc <- lfc[!duplicated(names(lfc))]
  lfc <- sort(lfc, decreasing = TRUE)
  
  assert_that(length(lfc) >= 200, "Ranked Entrez list is too small; mapping likely failed.")
  list(entrez_map = entrez_map, lfc_rank = lfc)
}

add_symbols <- function(entrez_ids) {
  AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = entrez_ids,
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  )
}

run_kegg_ora <- function(entrez_ids, organism = "hsa", pvalue_cutoff = 0.05) {
  clusterProfiler::enrichKEGG(
    gene = entrez_ids,
    organism = organism,
    keyType = "ncbi-geneid",
    pvalueCutoff = pvalue_cutoff
  )
}

run_go_gsea <- function(lfc_rank, ont = "BP") {
  clusterProfiler::gseGO(
    geneList = lfc_rank,
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    ont = ont,
    keyType = "ENTREZID",
    eps = 1e-300
  )
}

run_kegg_gsea <- function(lfc_rank, organism = "hsa") {
  clusterProfiler::gseKEGG(
    geneList = lfc_rank,
    organism = organism,
    keyType = "ncbi-geneid",
    pvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500
  )
}

go_emapplot_safe <- function(ego, showCategory = 10) {
  if (is_empty_enrich(ego)) return(NULL)
  ego_simplified <- clusterProfiler::simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
  ego_sim <- enrichplot::pairwise_termsim(ego_simplified)
  enrichplot::emapplot(ego_sim, showCategory = showCategory) + ggtitle("GO emapplot")
}

go_cnetplot_safe <- function(ego, showCategory = 5) {
  if (is_empty_enrich(ego)) return(NULL)
  enrichplot::cnetplot(ego, showCategory = showCategory) + ggtitle("GO cnetplot")
}
