# R/01_deseq2.R
suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
})

prefilter_dds <- function(dds, min_total = 1) {
  keep <- rowSums(DESeq2::counts(dds)) > min_total
  dds[keep, , drop = FALSE]
}

run_deseq2 <- function(dds) {
  DESeq2::DESeq(dds)
}

get_deseq_results <- function(dds, alpha = 0.05) {
  res <- DESeq2::results(dds, alpha = alpha)
  res <- res[order(res$pvalue), ]
  res
}

shrink_lfc_safe <- function(dds, coef, type = "apeglm") {
  out <- tryCatch(
    DESeq2::lfcShrink(dds, coef = coef, type = type),
    error = function(e) NULL
  )
  if (is.null(out)) DESeq2::results(dds) else out
}

make_transforms <- function(dds, blind = TRUE) {
  list(
    vsd = DESeq2::vst(dds, blind = blind),
    rld = DESeq2::rlog(dds, blind = blind)
  )
}
