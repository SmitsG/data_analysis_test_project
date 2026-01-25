# ==============================================
# Airway DESeq2 + GSEA Analysis
# ==============================================

# -----------------------------
# Load required libraries
# -----------------------------
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
library(enrichplot)
library(DOSE)
library(here)
library(vsn)
library(dplyr)

# -----------------------------
# Step 0: Settings & Parameters
# -----------------------------
pval_cutoff <- 0.05
fc_cutoff <- 1.5
top_n_genes <- 20
gsea_showCategory <- 15
ma_lfc_threshold <- 0.5
kegg_gsea_pvalue <- 0.2

# -----------------------------
# Step 1: Load required packages
# -----------------------------
required_packages <- c("airway","DESeq2","pheatmap","ggplot2","RColorBrewer",
                       "clusterProfiler","org.Hs.eg.db","SummarizedExperiment",
                       "EnhancedVolcano","KEGGREST","enrichplot","DOSE","here",
                       "vsn","dplyr")
for(pkg in required_packages) {
  if(!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# -----------------------------
# Step 2: Setup directories
# -----------------------------
dirs <- c("plots","results")
for(d in dirs){
  if(!dir.exists(here::here(d))) dir.create(here::here(d))
}

# -----------------------------
# Step 3: Load airway dataset
# -----------------------------
data("airway")
colData(airway)

# -----------------------------
# Step 4: Create DESeq2 dataset
# -----------------------------
dds <- DESeqDataSet(airway, design = ~ dex)
dds <- dds[rowSums(counts(dds)) > 1, ]

# -----------------------------
# Step 5: Run DESeq2 & results
# -----------------------------
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$pvalue), ]
write.csv(as.data.frame(resOrdered), here::here("results","Airway_DESeq2_results.csv"))

# -----------------------------
# Step 6: Transformations
# -----------------------------
vsd <- vst(dds, blind = TRUE)
rld <- rlog(dds, blind = TRUE)

# -----------------------------
# Step 7: PCA Plot Function
# -----------------------------
plot_pca <- function(rld_obj, group_var="dex", filename="PCA_plot.png"){
  pca_data <- plotPCA(rld_obj, intgroup=group_var, returnData=TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"))
  p <- ggplot(pca_data, aes(PC1, PC2, color = !!sym(group_var))) +
    geom_point(size=4) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle("PCA of Airway Samples") + theme_minimal()
  ggsave(here::here("plots", filename), plot=p, width=7, height=5)
  return(p)
}
pca_plot <- plot_pca(rld)

# -----------------------------
# Step 8: MA Plot Function
# -----------------------------
plot_ma <- function(res_obj, filename="MA_plot.png", main_title="MA Plot"){
  plotMA(res_obj, main=main_title)
  ggsave(here::here("plots", filename), width=7, height=5)
}
plot_ma(res)

# -----------------------------
# Step 9: Top Variable Genes Heatmap
# -----------------------------
topGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), top_n_genes)
mat <- assay(vsd)[topGenes, ]
mat <- mat - rowMeans(mat)
pheatmap(mat,
         annotation_col = as.data.frame(colData(dds)[, "dex", drop=FALSE]),
         main = "Top 20 Variable Genes Heatmap",
         color = colorRampPalette(rev(brewer.pal(9,"RdBu")))(255),
         filename = here::here("plots","Top20_VariableGenes_Heatmap.png"))

# -----------------------------
# Step 10: Volcano Plot
# -----------------------------
volcano_plot <- EnhancedVolcano(res,
                                lab = rownames(res),
                                x = 'log2FoldChange',
                                y = 'pvalue',
                                pCutoff = pval_cutoff,
                                FCcutoff = fc_cutoff,
                                title = "Volcano plot of DE genes")
ggsave(here::here("plots","Volcano_plot.png"), plot=volcano_plot, width=7, height=5)

# -----------------------------
# Step 11: Sample Distance Heatmap
# -----------------------------
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
pheatmap(sampleDistMatrix,
         col=colors,
         main="Sample-to-sample distance heatmap",
         filename = here::here("plots","Sample_distance_heatmap.png"))

# -----------------------------
# Step 12: Identify significant genes
# -----------------------------
sig_genes_up <- rownames(res)[which(res$padj < pval_cutoff & res$log2FoldChange > 0)]
sig_genes_down <- rownames(res)[which(res$padj < pval_cutoff & res$log2FoldChange < 0)]
write.csv(sig_genes_up, here::here("results","Upregulated_genes.csv"), row.names=FALSE)
write.csv(sig_genes_down, here::here("results","Downregulated_genes.csv"), row.names=FALSE)

# -----------------------------
# Step 13: GO Enrichment
# -----------------------------
ego <- enrichGO(gene = sig_genes_up,
                OrgDb = org.Hs.eg.db,
                keyType = "ENSEMBL",
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = pval_cutoff)
dotplot(ego, showCategory = 10) + ggtitle("GO Enrichment - Upregulated Genes")
ggsave(here::here("plots","GO_Enrichment_Dotplot.png"), width=7, height=5)

# -----------------------------
# Step 14: KEGG Enrichment
# -----------------------------
ensembl_genes_clean <- gsub("\\..*$","",rownames(res))
entrez_map <- bitr(ensembl_genes_clean, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
gene_list_entrez <- res$log2FoldChange[match(entrez_map$ENSEMBL, ensembl_genes_clean)]
names(gene_list_entrez) <- entrez_map$ENTREZID
gene_list_entrez <- sort(na.omit(gene_list_entrez[!duplicated(names(gene_list_entrez))]), decreasing=TRUE)

kegg_res <- enrichKEGG(gene = names(gene_list_entrez),
                       organism = "hsa", keyType="ncbi-geneid", pvalueCutoff=pval_cutoff)
dotplot(kegg_res, showCategory=10) + ggtitle("KEGG Pathway Enrichment")
ggsave(here::here("plots","KEGG_Enrichment_Dotplot.png"), width=7, height=5)

# -----------------------------
# Step 15: GSEA GO
# -----------------------------
gse_res <- gseGO(geneList = gene_list_entrez,
                 OrgDb = org.Hs.eg.db,
                 ont = "BP",
                 keyType = "ENTREZID",
                 eps=1e-300)
ridge_go <- ridgeplot(gse_res, showCategory=gsea_showCategory) +
  labs(x="Enrichment distribution", title="GO Enrichment Ridgeplot")
ggsave(here::here("plots","GSEA_GO_ridgeplot.png"), plot=ridge_go, width=8, height=6)

# -----------------------------
# Step 16: GSEA KEGG
# -----------------------------
gse_kegg <- gseKEGG(geneList = gene_list_entrez,
                    organism="hsa", keyType="ncbi-geneid",
                    pvalueCutoff=kegg_gsea_pvalue,
                    minGSSize=10, maxGSSize=500)
ridge_kegg <- ridgeplot(gse_kegg, showCategory=gsea_showCategory) +
  labs(x="Enrichment distribution", title="KEGG Enrichment Ridgeplot")
ggsave(here::here("plots","GSEA_KEGG_ridgeplot.png"), plot=ridge_kegg, width=8, height=6)

dot_go <- dotplot(gse_res, showCategory=10, split=".sign") + DOSE::facet_grid(.~.sign)
dot_kegg <- dotplot(gse_kegg, showCategory=10, split=".sign") + DOSE::facet_grid(.~.sign)
ggsave(here::here("plots","GSEA_GO_dotplot_sign.png"), plot=dot_go, width=8, height=6)
ggsave(here::here("plots","GSEA_KEGG_dotplot_sign.png"), plot=dot_kegg, width=8, height=6)

# -----------------------------
# Step 17: Mean-SD Plots
# -----------------------------
vsd_meanSD <- meanSdPlot(assay(vsd))$gg + ggtitle("Mean-SD Plot (VST)")
rld_meanSD <- meanSdPlot(assay(rld))$gg + ggtitle("Mean-SD Plot (RLD)")
ntd_meanSD <- meanSdPlot(assay(dds))$gg + ggtitle("Mean-SD Plot (untransformed)")
ggsave(here::here("plots","MeanSD_VST.png"), plot=vsd_meanSD, width=7, height=5)
ggsave(here::here("plots","MeanSD_RLD.png"), plot=rld_meanSD, width=7, height=5)
ggsave(here::here("plots","MeanSD_NTD.png"), plot=ntd_meanSD, width=7, height=5)

# -----------------------------
# Step 18: Cook's distance
# -----------------------------
png(here::here("plots","Cooks_distance_boxplot.png"), width=800, height=600)
boxplot(log10(assays(dds)[["cooks"]]), main="Log10 Cooks Distance")
dev.off()

# -----------------------------
# Step 19: Dispersion plot
# -----------------------------
png(here::here("plots","Dispersion_plot.png"), width=800, height=600)
plotDispEsts(dds, main="Dispersion Estimates")
dev.off()

# -----------------------------
# Step 20: Multi-panel MA plots
# -----------------------------
resGA <- results(dds, lfcThreshold=ma_lfc_threshold, altHypothesis="greaterAbs")
resLA <- results(dds, lfcThreshold=ma_lfc_threshold, altHypothesis="lessAbs")
resG  <- results(dds, lfcThreshold=ma_lfc_threshold, altHypothesis="greater")
resL  <- results(dds, lfcThreshold=ma_lfc_threshold, altHypothesis="less")
drawLines <- function() abline(h=c(-0.5,0.5), col="dodgerblue", lwd=2)
png(here::here("plots","MA_multi_panel.png"), width=1200, height=800)
par(mfrow=c(2,2), mar=c(2,2,1,1))
plotMA(resGA); drawLines()
plotMA(resLA); drawLines()
plotMA(resG); drawLines()
plotMA(resL); drawLines()
dev.off()

# -----------------------------
# Step 21: ORA full enrichment map
# -----------------------------
ego_full <- enrichGO(gene = rownames(res)[res$padj < pval_cutoff],
                     OrgDb=org.Hs.eg.db, keyType="ENSEMBL", ont="BP")
ego_sim_full <- pairwise_termsim(ego_full)
emap_full <- emapplot(ego_sim_full, showCategory=10) + ggtitle("Full ORA Enrichment Map")
ggsave(here::here("plots","ORA_full_enrichment_map.png"), plot=emap_full, width=10, height=8)

# -----------------------------
# Step 22: Save all key plots in single PDF
# -----------------------------
pdf(here::here("plots","Airway_DE_Analysis_All_Plots.pdf"), width=10, height=8)
print(volcano_plot)
pheatmap(mat,
         annotation_col = as.data.frame(colData(dds)[, "dex", drop=FALSE]),
         main="Top 20 Variable Genes Heatmap",
         color=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))
print(dotplot(ego, showCategory=10) + ggtitle("GO Enrichment - Upregulated Genes"))
print(dotplot(kegg_res, showCategory=10) + ggtitle("KEGG Pathway Enrichment"))
print(vsd_meanSD)
print(rld_meanSD)
print(ntd_meanSD)
print(ridge_go)
print(ridge_kegg)
print(dot_go)
print(dot_kegg)
print(emap_full)
dev.off()
