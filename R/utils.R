msg <- function(..., quiet = FALSE) {
  if (!isTRUE(quiet)) {
    message(...)
  }
}

require_namespace <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required R package is missing: ", pkg, call. = FALSE)
  }
}

parse_region <- function(region) {
  m <- regexec("^([^:]+):(\\d+)-(\\d+)$", region)
  x <- regmatches(region, m)[[1]]
  if (length(x) != 4) {
    stop("Region must look like chr5:12300000-12400000", call. = FALSE)
  }
  start <- as.integer(x[3])
  end <- as.integer(x[4])
  if (is.na(start) || is.na(end) || start > end) {
    stop("Invalid region coordinates: ", region, call. = FALSE)
  }
  list(chr = x[2], start = start, end = end)
}

sanitize_prefix <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

check_command <- function(cmd, label) {
  ok <- nzchar(Sys.which(cmd))
  if (!ok) {
    stop(label, " was not found in PATH: ", cmd, call. = FALSE)
  }
}

run_cmd <- function(command, args, quiet = FALSE) {
  if (!quiet) {
    message("$ ", command, " ", paste(args, collapse = " "))
  }
  status <- system2(
    command,
    args = args,
    stdout = if (quiet) FALSE else "",
    stderr = if (quiet) FALSE else ""
  )
  if (!identical(status, 0L)) {
    stop("Command failed with status ", status, ": ", command, call. = FALSE)
  }
  invisible(TRUE)
}

find_beagle <- function(beagle_jar = NULL) {
  candidates <- c(
    beagle_jar,
    Sys.getenv("BEAGLE_JAR", unset = NA_character_),
    file.path("tools", "beagle.jar")
  )
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) normalizePath(hit[1], winslash = "/", mustWork = TRUE) else NULL
}

format_beagle_args <- function(x) {
  if (!length(x)) {
    return(character())
  }
  if (is.null(names(x)) || any(!nzchar(names(x)))) {
    stop("beagle_args must be a named list, for example list(impute = FALSE)",
         call. = FALSE)
  }
  vapply(names(x), function(nm) {
    val <- x[[nm]]
    if (is.logical(val)) {
      val <- tolower(as.character(val))
    }
    paste0(nm, "=", val)
  }, character(1))
}

open_text <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
}

read_vcf_genotypes <- function(path,
                               region_info = NULL,
                               snp_only = TRUE,
                               allow_unphased = FALSE) {
  con <- open_text(path)
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  header_i <- grep("^#CHROM\\t", lines)
  if (!length(header_i)) {
    stop("VCF header line was not found: ", path, call. = FALSE)
  }
  header <- strsplit(lines[header_i[1]], "\t", fixed = TRUE)[[1]]
  if (length(header) < 10) {
    stop("VCF has no sample genotype columns: ", path, call. = FALSE)
  }
  samples <- header[10:length(header)]
  var_lines <- lines[!grepl("^#", lines)]
  if (!length(var_lines)) {
    stop("No variants found in regional VCF: ", path, call. = FALSE)
  }

  records <- strsplit(var_lines, "\t", fixed = TRUE)
  keep <- rep(TRUE, length(records))
  if (!is.null(region_info)) {
    keep <- vapply(records, function(z) {
      z[1] == region_info$chr &&
        suppressWarnings(as.integer(z[2])) >= region_info$start &&
        suppressWarnings(as.integer(z[2])) <= region_info$end
    }, logical(1))
  }
  records <- records[keep]
  if (!length(records)) {
    stop("No variants remain after applying region filter.", call. = FALSE)
  }

  if (snp_only) {
    is_snp <- vapply(records, function(z) {
      ref <- z[4]
      alt <- strsplit(z[5], ",", fixed = TRUE)[[1]]
      nchar(ref) == 1 && all(nchar(alt) == 1) &&
        grepl("^[ACGTNacgtn]$", ref) &&
        all(grepl("^[ACGTNacgtn]$", alt))
    }, logical(1))
    records <- records[is_snp]
  }
  if (!length(records)) {
    stop("No SNP variants remain for haplotype network construction.",
         call. = FALSE)
  }

  nvar <- length(records)
  hap1 <- matrix("N", nrow = length(samples), ncol = nvar,
                 dimnames = list(samples, NULL))
  hap2 <- matrix("N", nrow = length(samples), ncol = nvar,
                 dimnames = list(samples, NULL))
  variants <- data.frame(
    marker = character(nvar),
    chrom = character(nvar),
    pos = integer(nvar),
    id = character(nvar),
    ref = character(nvar),
    alt = character(nvar),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(records)) {
    z <- records[[i]]
    fmt <- strsplit(z[9], ":", fixed = TRUE)[[1]]
    gt_idx <- match("GT", fmt)
    if (is.na(gt_idx)) {
      stop("FORMAT column lacks GT at variant ", z[1], ":", z[2],
           call. = FALSE)
    }
    alleles <- toupper(c(z[4], strsplit(z[5], ",", fixed = TRUE)[[1]]))
    variants[i, ] <- list(
      marker = paste0(z[1], ":", z[2]),
      chrom = z[1],
      pos = as.integer(z[2]),
      id = z[3],
      ref = z[4],
      alt = z[5]
    )
    for (j in seq_along(samples)) {
      gt_field <- strsplit(z[9 + j], ":", fixed = TRUE)[[1]]
      gt <- gt_field[gt_idx]
      parsed <- parse_gt(gt, alleles, allow_unphased = allow_unphased)
      hap1[j, i] <- parsed[1]
      hap2[j, i] <- parsed[2]
    }
  }

  colnames(hap1) <- variants$marker
  colnames(hap2) <- variants$marker
  list(samples = samples, variants = variants, hap1 = hap1, hap2 = hap2)
}

parse_gt <- function(gt, alleles, allow_unphased = FALSE) {
  if (is.na(gt) || gt %in% c(".", "./.", ".|.")) {
    return(c("N", "N"))
  }
  if (grepl("/", gt, fixed = TRUE) && !allow_unphased) {
    stop(
      "Unphased genotype found (", gt, "). Run Beagle or set ",
      "allow_unphased = TRUE for already ordered demo data.",
      call. = FALSE
    )
  }
  parts <- strsplit(gt, "[|/]")[[1]]
  if (length(parts) != 2) {
    stop("Only diploid genotypes are supported; got GT=", gt, call. = FALSE)
  }
  vapply(parts, function(a) {
    if (a == ".") {
      return("N")
    }
    idx <- suppressWarnings(as.integer(a)) + 1L
    if (is.na(idx) || idx < 1L || idx > length(alleles)) {
      return("N")
    }
    base <- alleles[idx]
    if (!grepl("^[ACGTN]$", base)) "N" else base
  }, character(1))
}

read_metadata <- function(metadata, samples, sample_col = NULL, group_col = NULL) {
  if (is.null(metadata) || !nzchar(metadata)) {
    return(data.frame(sample = samples, group = "All", stringsAsFactors = FALSE))
  }
  if (!file.exists(metadata)) {
    stop("Metadata file does not exist: ", metadata, call. = FALSE)
  }
  first <- readLines(metadata, n = 1, warn = FALSE)
  sep <- if (grepl(",", first, fixed = TRUE)) "," else ""
  meta <- utils::read.table(
    metadata,
    header = TRUE,
    sep = sep,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = "",
    quote = "\""
  )
  if (!nrow(meta)) {
    stop("Metadata has no rows: ", metadata, call. = FALSE)
  }
  if (is.null(sample_col)) {
    candidates <- c("sample", "sample_id", "Sample", "SampleID", "id", "ID")
    sample_col <- candidates[candidates %in% names(meta)][1]
    if (is.na(sample_col)) {
      sample_col <- names(meta)[1]
    }
  }
  if (!sample_col %in% names(meta)) {
    stop("sample_col not found in metadata: ", sample_col, call. = FALSE)
  }
  if (is.null(group_col)) {
    rest <- setdiff(names(meta), sample_col)
    if (!length(rest)) {
      group_col <- sample_col
    } else {
      group_col <- rest[1]
    }
  }
  if (!group_col %in% names(meta)) {
    stop("group_col not found in metadata: ", group_col, call. = FALSE)
  }
  out <- data.frame(
    sample = as.character(meta[[sample_col]]),
    group = as.character(meta[[group_col]]),
    stringsAsFactors = FALSE
  )
  out$group[is.na(out$group) | !nzchar(out$group)] <- "Unknown"
  out <- out[!duplicated(out$sample), , drop = FALSE]
  missing <- setdiff(samples, out$sample)
  if (length(missing)) {
    out <- rbind(
      out,
      data.frame(sample = missing, group = "Unknown", stringsAsFactors = FALSE)
    )
  }
  out[match(samples, out$sample), , drop = FALSE]
}

build_haplotypes <- function(vcf_data, metadata) {
  samples <- vcf_data$samples
  hap_mat <- rbind(vcf_data$hap1, vcf_data$hap2)
  rownames(hap_mat) <- c(paste0(samples, "_1"), paste0(samples, "_2"))
  group <- rep(metadata$group, 2)
  sample <- rep(samples, 2)
  copy <- rep(c("hap1", "hap2"), each = length(samples))

  dna <- ape::as.DNAbin(hap_mat)
  hap <- pegas::haplotype(dna)
  hap_labels <- paste0("H", formatC(seq_len(nrow(hap)), width = 2, flag = "0"))
  rownames(hap) <- hap_labels
  index <- attr(hap, "index")

  sample_haplotypes <- data.frame(
    sample = sample,
    copy = copy,
    group = group,
    hap_id = NA_character_,
    stringsAsFactors = FALSE
  )
  for (i in seq_along(index)) {
    sample_haplotypes$hap_id[index[[i]]] <- hap_labels[i]
  }
  sample_haplotypes$sequence <- apply(
    hap_mat[rownames(hap_mat) %in% rownames(hap_mat), , drop = FALSE],
    1,
    paste0,
    collapse = ""
  )

  counts <- table(factor(sample_haplotypes$hap_id, levels = hap_labels))
  haplotypes <- data.frame(
    hap_id = hap_labels,
    count = as.integer(counts),
    frequency = as.integer(counts) / sum(counts),
    sequence = vapply(index, function(i) {
      paste0(as.character(hap_mat[i[1], ]), collapse = "")
    }, character(1)),
    stringsAsFactors = FALSE
  )

  frequency <- as.data.frame(
    table(
      hap_id = factor(sample_haplotypes$hap_id, levels = hap_labels),
      group = sample_haplotypes$group
    ),
    stringsAsFactors = FALSE
  )
  names(frequency)[3] <- "count"
  group_totals <- tapply(frequency$count, frequency$group, sum)
  frequency$within_group_frequency <- frequency$count /
    group_totals[frequency$group]
  frequency$total_frequency <- frequency$count / sum(frequency$count)
  frequency <- frequency[frequency$count > 0, , drop = FALSE]

  group_levels <- sort(unique(sample_haplotypes$group))
  pie <- matrix(0, nrow = length(hap_labels), ncol = length(group_levels),
                dimnames = list(hap_labels, group_levels))
  for (i in seq_along(index)) {
    pie[i, ] <- table(factor(group[index[[i]]], levels = group_levels))
  }

  net <- pegas::haploNet(hap)
  list(
    hap = hap,
    net = net,
    pie = pie,
    haplotypes = haplotypes,
    frequency = frequency,
    sample_haplotypes = sample_haplotypes
  )
}

draw_hapnet <- function(hap,
                        net,
                        pie,
                        freq,
                        region,
                        outdir,
                        prefix,
                        export,
                        label,
                        min_hap_count,
                        palette,
                        width,
                        height,
                        dpi) {
  pal <- get_palette(palette, ncol(pie))
  names(pal) <- colnames(pie)
  node_sizes <- scale_size(sqrt(freq), to = c(0.8, 5.5))
  label_flag <- isTRUE(label)
  if (identical(label, "auto")) {
    label_flag <- sum(freq >= min_hap_count) <= 30
  }

  draw_one <- function() {
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op), add = TRUE)
    graphics::par(
      mar = c(0.5, 0.5, 1.8, 0.5),
      bg = "white",
      fg = "#202020",
      family = "sans",
      xpd = NA
    )
    plot(
      net,
      size = node_sizes,
      pie = pie,
      bg = pal,
      col = "#202020",
      lwd = 1.25,
      labels = label_flag,
      cex = 0.82,
      font = 2,
      show.mutation = 2,
      threshold = 0,
      fast = TRUE
    )
    graphics::title(main = paste0(region, " haplotype network"), cex.main = 0.95,
                    font.main = 2, line = 0.5)
    graphics::legend(
      "bottomleft",
      legend = colnames(pie),
      fill = pal,
      border = "#202020",
      bty = "n",
      cex = 0.82,
      inset = 0.01
    )
  }

  files <- character()
  if ("pdf" %in% export) {
    f <- file.path(outdir, paste0(prefix, ".hapnet.pdf"))
    grDevices::cairo_pdf(f, width = width, height = height, onefile = TRUE)
    draw_one()
    grDevices::dev.off()
    files["pdf"] <- normalizePath(f, winslash = "/", mustWork = FALSE)
  }
  if ("svg" %in% export) {
    f <- file.path(outdir, paste0(prefix, ".hapnet.svg"))
    if (requireNamespace("svglite", quietly = TRUE)) {
      svglite::svglite(f, width = width, height = height)
    } else {
      grDevices::svg(f, width = width, height = height)
    }
    draw_one()
    grDevices::dev.off()
    files["svg"] <- normalizePath(f, winslash = "/", mustWork = FALSE)
  }
  if ("tiff" %in% export) {
    f <- file.path(outdir, paste0(prefix, ".hapnet.tiff"))
    grDevices::tiff(
      f,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      compression = "lzw",
      type = "cairo"
    )
    draw_one()
    grDevices::dev.off()
    files["tiff"] <- normalizePath(f, winslash = "/", mustWork = FALSE)
  }
  files
}

get_palette <- function(palette, n) {
  base <- switch(
    palette,
    nature = c(
      "#3B4992", "#EE0000", "#008B45", "#631879", "#008280",
      "#BB0021", "#5F559B", "#A20056", "#808180", "#1B1919"
    ),
    okabe_ito = c(
      "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
      "#56B4E9", "#F0E442", "#000000"
    ),
    stop("Unknown palette: ", palette, call. = FALSE)
  )
  if (n <= length(base)) {
    base[seq_len(n)]
  } else {
    grDevices::colorRampPalette(base)(n)
  }
}

scale_size <- function(x, to = c(0.8, 5.5)) {
  if (!length(x)) {
    return(x)
  }
  if (max(x) == min(x)) {
    return(rep(mean(to), length(x)))
  }
  to[1] + (x - min(x)) / (max(x) - min(x)) * diff(to)
}

write_tables <- function(built, variants, outdir, prefix) {
  files <- c(
    haplotypes = file.path(outdir, paste0(prefix, ".haplotypes.tsv")),
    frequency = file.path(outdir, paste0(prefix, ".haplotype_frequency.tsv")),
    counts_wide = file.path(outdir, paste0(prefix, ".haplotype_counts_wide.tsv")),
    sample_haplotypes = file.path(outdir, paste0(prefix, ".sample_haplotypes.tsv")),
    variants = file.path(outdir, paste0(prefix, ".variants.tsv"))
  )
  utils::write.table(
    built$haplotypes, files["haplotypes"], sep = "\t", quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    built$frequency, files["frequency"], sep = "\t", quote = FALSE,
    row.names = FALSE
  )
  wide <- as.data.frame.matrix(built$pie)
  wide <- cbind(hap_id = rownames(wide), wide)
  utils::write.table(
    wide, files["counts_wide"], sep = "\t", quote = FALSE, row.names = FALSE
  )
  utils::write.table(
    built$sample_haplotypes, files["sample_haplotypes"], sep = "\t",
    quote = FALSE, row.names = FALSE
  )
  utils::write.table(
    variants, files["variants"], sep = "\t", quote = FALSE, row.names = FALSE
  )
  normalizePath(files, winslash = "/", mustWork = FALSE)
}
