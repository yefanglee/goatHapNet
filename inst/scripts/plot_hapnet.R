#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!grepl("^--", arg)) {
      next
    }
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- kv[1]
    val <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else TRUE
    out[[key]] <- val
  }
  out
}

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0 || is.na(a) || !nzchar(a)) b else a
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
need <- c("vcf", "region")
missing <- need[!need %in% names(args)]
if (length(missing)) {
  stop("Missing required arguments: ", paste(missing, collapse = ", "),
       call. = FALSE)
}

library(goatHapNet)

as_logical <- function(x, default = FALSE) {
  if (is.null(x)) {
    return(default)
  }
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

res <- plot_hapnet(
  vcf = args$vcf,
  region = args$region,
  metadata = args$metadata %||% NULL,
  outdir = args$outdir %||% "hapnet_output",
  prefix = args$prefix %||% NULL,
  sample_col = args$sample_col %||% NULL,
  group_col = args$group_col %||% NULL,
  extract_region = as_logical(args$extract_region, default = TRUE),
  run_beagle = as_logical(args$run_beagle, default = TRUE),
  beagle_jar = args$beagle_jar %||% NULL,
  bcftools = args$bcftools %||% "bcftools",
  java = args$java %||% "java",
  threads = as.integer(args$threads %||% 4),
  java_mem = args$java_mem %||% "8g",
  allow_unphased = as_logical(args$allow_unphased, default = FALSE)
)

cat("Wrote files:\n")
cat(paste0("  ", unname(res$files), collapse = "\n"), "\n")
