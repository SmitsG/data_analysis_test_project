# _targets.R
library(targets)
library(tarchetypes)
library(here)

# Keep objects out of global env
tar_option_set(
  packages = c(
    "airway","DESeq2","SummarizedExperiment","ggplot2","pheatmap","RColorBrewer",
    "clusterProfiler","org.Hs.eg.db","enrichplot","DOSE","EnhancedVolcano","vsn",
    "dplyr","matrixStats","AnnotationDbi","here"
  ),
  format = "rds"
)

source(here::here("R", "functions.R"))

list(
  tar_target(dirs, ensure_dirs(), format = "file"),
  
  # Data
  tar_target(airway_se, { data("airway", package = "airway"); airway }),
  
  # DESeq2
  tar_target(dds, make_dds(airway_se)),
  tar_target(dds_fit, run_deseq(dds)),
  tar_target(res_list, get_results(dds_fit)),
  tar_target(res, res_list$res),
  tar_target(resOrdered, res_list$resOrdered),
  
  # Save results
  tar_target(res_csv, save_csv(as.data.frame(resOrdered), here::here("results","Airway_DESeq2_results.csv")), format = "file"),
  
  # Transformations
  tar_target(tx, make_transforms(dds_fit)),
  tar_target(vsd, tx$vsd),
  tar_target(rld, tx$rld),
  
  # PCA
  tar_target(pca_plot, plot_pca(rld)),
  tar_target(pca_png, save_gg(pca_plot, here::here("plots","PCA_plot.png")), format = "file"),
  
  # MA plot
  tar_target(ma_png, plot_ma_png(res, here::here("plots","MA_plot.png")), format = "file"),
  
  # Heatmaps
  tar_target(topvar_heatmap_png, heatmap_top_var(vsd, dds_fit, n = 20, path = here::here("plots","Top20_VariableGenes_Heatmap.png")), format = "file"),
  tar_target(sampledist_heatmap_png, sample_distance_heatmap(vsd, here::here("plots","Sample_distance_heatmap.png")), format = "file"),
  
  # Volcano
  tar_target(volcano_plot, plot_volcano(res)),
  tar_target(volcano_png, save_gg(volcano_plot, here::here("plots","Volcano_plot.png")), format = "file"),
  
  # Significant genes
  tar_target(sig, get_sig_genes(res)),
  tar_target(up_csv, save_csv(data.frame(ENSEMBL = sig$up), here::here("results","Upregulated_genes.csv")), format = "file"),
  tar_target(down_csv, save_csv(data.frame(ENSEMBL = sig$down), here::here("results","Downregulated_genes.csv")), format = "file"),
  
  # GO ORA
  tar_target(ego, ora_go_up(sig$up)),
  tar_target(go_dot, dotplot(ego, showCategory = 10) + ggtitle("GO Enrichment - Upregulated Genes")),
  tar_target(go_dot_png, save_gg(go_dot, here::here("plots","GO_Enrichment_Dotplot.png")), format = "file"),
  
  # Mapping + KEGG ORA
  tar_target(map_list, map_ensembl_to_entrez(res)),
  tar_target(gene_list_entrez, map_list$gene_list_entrez),
  tar_target(kegg_res, ora_kegg(gene_list_entrez)),
  tar_target(kegg_dot, dotplot(kegg_res, showCategory = 10) + ggtitle("KEGG Pathway Enrichment")),
  tar_target(kegg_dot_png, save_gg(kegg_dot, here::here("plots","KEGG_Enrichment_Dotplot.png")), format = "file"),
  
  # Save ranked gene lists (basic + symbols)
  tar_target(
    ranked_gene_csvs,
    save_ranked_genes(
      gene_list_entrez,
      here::here("results","Ranked_genes_for_GSEA.csv"),
      here::here("results","Ranked_Genes_for_GSEA_with_SYMBOL.csv")
    ),
    format = "file"
  ),
  
  # GSEA
  tar_target(gse_go, gsea_go(gene_list_entrez)),
  tar_target(gse_kegg, gsea_kegg(gene_list_entrez)),
  
  tar_target(ridge_go, ridgeplot(gse_go, showCategory = 15) + labs(x="Enrichment distribution", title="GO Enrichment Ridgeplot")),
  tar_target(ridge_kegg, ridgeplot(gse_kegg, showCategory = 15) + labs(x="Enrichment distribution", title="KEGG Enrichment Ridgeplot")),
  tar_target(ridge_go_png, save_gg(ridge_go, here::here("plots","GSEA_GO_ridgeplot.png"), width = 8, height = 6), format = "file"),
  tar_target(ridge_kegg_png, save_gg(ridge_kegg, here::here("plots","GSEA_KEGG_ridgeplot.png"), width = 8, height = 6), format = "file"),
  
  tar_target(dot_go, dotplot(gse_go, showCategory = 10, split = ".sign") + DOSE::facet_grid(.~.sign)),
  tar_target(dot_kegg, dotplot(gse_kegg, showCategory = 10, split = ".sign") + DOSE::facet_grid(.~.sign)),
  tar_target(dot_go_png, save_gg(dot_go, here::here("plots","GSEA_GO_dotplot_sign.png"), width = 8, height = 6), format = "file"),
  tar_target(dot_kegg_png, save_gg(dot_kegg, here::here("plots","GSEA_KEGG_dotplot_sign.png"), width = 8, height = 6), format = "file"),
  
  # Extra QC
  tar_target(mean_sd, mean_sd_plots(vsd, rld, dds_fit)),
  tar_target(meansd_vst_png, save_gg(mean_sd$vsd, here::here("plots","MeanSD_VST.png")), format = "file"),
  tar_target(meansd_rld_png, save_gg(mean_sd$rld, here::here("plots","MeanSD_RLD.png")), format = "file"),
  tar_target(meansd_ntd_png, save_gg(mean_sd$ntd, here::here("plots","MeanSD_NTD.png")), format = "file"),
  tar_target(cooks_png, cooks_png(dds_fit, here::here("plots","Cooks_distance_boxplot.png")), format = "file"),
  tar_target(disp_png, dispersion_png(dds_fit, here::here("plots","Dispersion_plot.png")), format = "file"),
  tar_target(ma_multi_png, ma_multipanel_png(dds_fit, here::here("plots","MA_multi_panel.png")), format = "file"),
  
  # Enrichment network / maps
  tar_target(go_cnet, go_cnetplot(ego, gene_list_entrez)),
  tar_target(go_cnet_png, save_gg(go_cnet, here::here("plots","GO_cnetplot.png"), width = 10, height = 8), format = "file"),
  
  tar_target(go_emap, go_emapplot(ego)),
  tar_target(go_emap_png, save_gg(go_emap, here::here("plots","GO_emapplot.png"), width = 10, height = 8), format = "file"),
  
  tar_target(ora_full_map, ora_full_emap(res)),
  tar_target(ora_full_map_png, save_gg(ora_full_map, here::here("plots","ORA_full_enrichment_map.png"), width = 10, height = 8), format = "file"),
  
  # PDF report (single artifact)
  tar_target(
    pdf_report,
    save_pdf_report(
      here::here("plots","Airway_DE_Analysis_All_Plots.pdf"),
      res = res,
      dds = dds_fit,
      vsd = vsd,
      ego = ego,
      kegg_res = kegg_res,
      mean_sd = mean_sd,
      ridge_go = ridge_go,
      ridge_kegg = ridge_kegg,
      dot_go = dot_go,
      dot_kegg = dot_kegg,
      emap_full = ora_full_map
    ),
    format = "file"
  )
)
