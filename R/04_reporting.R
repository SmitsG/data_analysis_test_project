# R/04_reporting.R
suppressPackageStartupMessages({
  library(ggplot2)
  library(enrichplot)
  library(DOSE)
})

safe_dotplot <- function(enrich_obj, title, showCategory = 10) {
  if (is_empty_enrich(enrich_obj)) return(NULL)
  enrichplot::dotplot(enrich_obj, showCategory = showCategory) + ggtitle(title)
}

safe_ridgeplot <- function(gsea_obj, title, showCategory = 15) {
  if (is_empty_enrich(gsea_obj)) return(NULL)
  enrichplot::ridgeplot(gsea_obj, showCategory = showCategory) + ggtitle(title)
}

safe_dotplot_split <- function(gsea_obj, title, showCategory = 10) {
  if (is_empty_enrich(gsea_obj)) return(NULL)
  enrichplot::dotplot(gsea_obj, showCategory = showCategory, split = ".sign") +
    DOSE::facet_grid(. ~ .sign) +
    ggtitle(title)
}

write_pdf_report <- function(pdf_path, ggplots = list(), base_plot_fns = list()) {
  grDevices::pdf(pdf_path, width = 10, height = 8)
  on.exit(grDevices::dev.off(), add = TRUE)
  
  for (p in ggplots) if (!is.null(p)) print(p)
  for (f in base_plot_fns) if (is.function(f)) f()
  
  pdf_path
}
