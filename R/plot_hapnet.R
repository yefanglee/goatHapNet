#' Plot a haplotype network from a VCF region
#'
#' `plot_hapnet()` extracts a genomic interval, optionally phases it with
#' Beagle, builds haplotypes, computes group frequencies, and writes
#' publication-ready PDF/SVG/TIFF haplotype network plots.
#'
#' @param vcf Input VCF path. Gzipped VCF is recommended.
#' @param region Region string such as `"chr5:12300000-12400000"`.
#' @param metadata Optional sample metadata file. Tab, comma, or whitespace
#'   delimited files are accepted.
#' @param outdir Output directory.
#' @param prefix Output prefix. If `NULL`, a prefix is derived from `region`.
#' @param sample_col Metadata sample ID column. If `NULL`, inferred.
#' @param group_col Metadata group column used for pie charts. If `NULL`,
#'   inferred from the first non-sample column.
#' @param extract_region Extract `region` from `vcf` with bcftools. Set to
#'   `FALSE` if `vcf` is already a regional VCF.
#' @param run_beagle Run Beagle phasing after extracting the region.
#' @param beagle_jar Path to Beagle jar. If `NULL`, checks `BEAGLE_JAR` and
#'   `tools/beagle.jar`.
#' @param beagle_args Named list of additional Beagle arguments.
#' @param bcftools `bcftools` executable.
#' @param java `java` executable.
#' @param threads Number of threads for Beagle.
#' @param java_mem Java heap, for example `"8g"`.
#' @param snp_only Keep only SNP variants for network construction.
#' @param allow_unphased Allow unphased genotypes. Recommended only for demo
#'   data or known phased inputs using `/` separators.
#' @param min_hap_count Minimum total haplotype-copy count considered when
#'   deciding whether `label = "auto"` should label the network.
#' @param label `TRUE`, `FALSE`, or `"auto"`.
#' @param palette Named palette. Currently `"nature"` or `"okabe_ito"`.
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#' @param dpi TIFF resolution.
#' @param export Output formats: any of `"pdf"`, `"svg"`, `"tiff"`.
#' @param keep_intermediate Keep intermediate VCF files.
#' @param quiet Suppress progress messages.
#'
#' @return Invisibly returns a list containing output files, haplotype tables,
#'   the `pegas` haplotype object, and the network object.
#' @export
plot_hapnet <- function(vcf,
                        region,
                        metadata = NULL,
                        outdir = "hapnet_output",
                        prefix = NULL,
                        sample_col = NULL,
                        group_col = NULL,
                        extract_region = TRUE,
                        run_beagle = TRUE,
                        beagle_jar = NULL,
                        beagle_args = list(),
                        bcftools = "bcftools",
                        java = "java",
                        threads = 4,
                        java_mem = "8g",
                        snp_only = TRUE,
                        allow_unphased = FALSE,
                        min_hap_count = 1,
                        label = "auto",
                        palette = "nature",
                        width = 7,
                        height = 6,
                        dpi = 600,
                        export = c("pdf", "svg", "tiff"),
                        keep_intermediate = TRUE,
                        quiet = FALSE) {
  require_namespace("ape")
  require_namespace("pegas")

  if (!file.exists(vcf)) {
    stop("VCF does not exist: ", vcf, call. = FALSE)
  }
  region_info <- parse_region(region)
  if (is.null(prefix)) {
    prefix <- sanitize_prefix(region)
  }
  export <- match.arg(export, c("pdf", "svg", "tiff"), several.ok = TRUE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  intdir <- file.path(outdir, "intermediate")
  dir.create(intdir, recursive = TRUE, showWarnings = FALSE)

  regional_vcf <- file.path(intdir, paste0(prefix, ".region.vcf.gz"))
  phased_vcf <- regional_vcf

  if (extract_region) {
    check_command(bcftools, "bcftools")
    msg("Extracting region with bcftools: ", region, quiet = quiet)
    run_cmd(
      bcftools,
      c("view", "-r", region, "-Oz", "-o", regional_vcf, vcf),
      quiet = quiet
    )
    run_cmd(bcftools, c("index", "-t", regional_vcf), quiet = quiet)
  } else {
    regional_vcf <- vcf
    phased_vcf <- regional_vcf
    msg("Skipping bcftools extraction because extract_region = FALSE",
        quiet = quiet)
  }

  if (run_beagle) {
    check_command(java, "java")
    beagle_jar <- find_beagle(beagle_jar)
    if (is.null(beagle_jar)) {
      stop(
        "Beagle jar was not found. Set BEAGLE_JAR, place tools/beagle.jar, ",
        "or pass beagle_jar = '/path/to/beagle.jar'.",
        call. = FALSE
      )
    }
    msg("Phasing regional VCF with Beagle", quiet = quiet)
    beagle_out_prefix <- file.path(intdir, paste0(prefix, ".beagle"))
    args <- c(
      paste0("-Xmx", java_mem),
      "-jar", beagle_jar,
      paste0("gt=", regional_vcf),
      paste0("out=", beagle_out_prefix),
      paste0("nthreads=", threads),
      format_beagle_args(beagle_args)
    )
    run_cmd(java, args, quiet = quiet)
    phased_vcf <- paste0(beagle_out_prefix, ".vcf.gz")
    if (!file.exists(phased_vcf)) {
      stop("Beagle finished but phased VCF was not found: ", phased_vcf,
           call. = FALSE)
    }
    if (extract_region) {
      run_cmd(bcftools, c("index", "-t", phased_vcf), quiet = quiet)
    }
  } else {
    allow_unphased <- isTRUE(allow_unphased)
  }

  msg("Parsing phased genotypes and building haplotypes", quiet = quiet)
  vcf_data <- read_vcf_genotypes(
    phased_vcf,
    region_info = region_info,
    snp_only = snp_only,
    allow_unphased = allow_unphased
  )
  meta <- read_metadata(metadata, vcf_data$samples, sample_col, group_col)
  built <- build_haplotypes(vcf_data, meta)

  msg("Drawing haplotype network", quiet = quiet)
  plot_files <- draw_hapnet(
    hap = built$hap,
    net = built$net,
    pie = built$pie,
    freq = built$haplotypes$count,
    region = region,
    outdir = outdir,
    prefix = prefix,
    export = export,
    label = label,
    min_hap_count = min_hap_count,
    palette = palette,
    width = width,
    height = height,
    dpi = dpi
  )

  table_files <- write_tables(
    built = built,
    variants = vcf_data$variants,
    outdir = outdir,
    prefix = prefix
  )

  if (!keep_intermediate) {
    unlink(intdir, recursive = TRUE, force = TRUE)
  }

  result <- list(
    files = c(plot_files, table_files),
    phased_vcf = if (keep_intermediate) phased_vcf else NA_character_,
    haplotypes = built$haplotypes,
    frequency = built$frequency,
    sample_haplotypes = built$sample_haplotypes,
    variants = vcf_data$variants,
    hap = built$hap,
    network = built$net
  )
  class(result) <- "goat_hapnet_result"
  invisible(result)
}
