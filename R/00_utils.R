# R/00_utils.R
suppressPackageStartupMessages({
  library(here)
})

init_output_dirs <- function(base_dir = here::here(),
                             plot_dir = "plots",
                             results_dir = "results") {
  stopifnot(is.character(base_dir), length(base_dir) == 1)
  
  plots_path <- file.path(base_dir, plot_dir)
  results_path <- file.path(base_dir, results_dir)
  
  dir.create(plots_path, showWarnings = FALSE, recursive = TRUE)
  dir.create(results_path, showWarnings = FALSE, recursive = TRUE)
  
  list(plots = plots_path, results = results_path)
}

safe_remove <- function(path) {
  if (file.exists(path)) file.remove(path)
  invisible(TRUE)
}

log_line <- function(path, ...) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  txt <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = " "))
  cat(txt, "\n")
  cat(txt, "\n", file = path, append = TRUE)
  invisible(TRUE)
}

assert_that <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
  invisible(TRUE)
}

save_csv <- function(df, path) {
  utils::write.csv(df, file = path, row.names = FALSE)
  path
}

save_gg <- function(p, path, width = 7, height = 5) {
  ggplot2::ggsave(filename = path, plot = p, width = width, height = height)
  path
}

save_base_png <- function(path, expr, width = 1000, height = 800, res = 120) {
  grDevices::png(path, width = width, height = height, res = res)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
  path
}

is_empty_enrich <- function(x) {
  is.null(x) || nrow(as.data.frame(x)) == 0
}
