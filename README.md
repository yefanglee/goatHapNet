# goatHapNet

`goatHapNet` is a small R package for automated haplotype network plotting from
candidate gene regions in goat VCF data.

The main entry point is:

```r
plot_hapnet(
  vcf = "example.vcf.gz",
  region = "chr5:12300000-12400000",
  metadata = "metadata.txt"
)
```

It performs the full workflow:

1. Extract the requested interval with `bcftools view`.
2. Phase the regional VCF with Beagle.
3. Parse phased genotypes into diploid haplotype copies.
4. Calculate haplotype frequencies per sample group.
5. Draw a publication-style haplotype network with pie charts.
6. Export `PDF`, `SVG`, and high-resolution `TIFF`.
7. Save intermediate tables for reproducibility.

## Requirements

R packages:

```r
install.packages(c("ape", "pegas", "svglite", "testthat"))
```

External tools:

- `bcftools` in `PATH`
- `java` in `PATH`
- Beagle jar, provided either by:
  - setting `BEAGLE_JAR=/path/to/beagle.jar`, or
  - placing it at `tools/beagle.jar`, or
  - passing `beagle_jar = "/path/to/beagle.jar"` to `plot_hapnet()`

## Install

From the parent directory:

```bash
R CMD INSTALL goatHapNet
```

Or during development:

```r
remotes::install_local(".")
```

## Metadata Format

The metadata file should be tab-delimited or comma-delimited. The first column is
used as the sample ID unless `sample_col` is provided. The first non-sample
column is used as the grouping variable unless `group_col` is provided.

Example:

```text
sample  population
G001    North
G002    North
G003    South
G004    South
```

Sample names must match VCF sample names.

## Real Data Example

```r
library(goatHapNet)

res <- plot_hapnet(
  vcf = "goats.allchr.vcf.gz",
  region = "chr5:12300000-12400000",
  metadata = "metadata.txt",
  group_col = "population",
  outdir = "hapnet_chr5_geneA",
  prefix = "geneA_chr5",
  threads = 8,
  export = c("pdf", "svg", "tiff")
)

res$files
res$haplotypes
res$frequency
```

## Demo Without Beagle

The bundled example VCF is already phased, so it can be used without Beagle:

```r
library(goatHapNet)

plot_hapnet(
  vcf = system.file("extdata", "example.vcf", package = "goatHapNet"),
  region = "chr5:12300000-12400000",
  metadata = system.file("extdata", "metadata.txt", package = "goatHapNet"),
  group_col = "population",
  extract_region = FALSE,
  run_beagle = FALSE,
  outdir = "demo_hapnet"
)
```

## Outputs

For a prefix such as `chr5_12300000_12400000`, the output directory contains:

- `chr5_12300000_12400000.hapnet.pdf`
- `chr5_12300000_12400000.hapnet.svg`
- `chr5_12300000_12400000.hapnet.tiff`
- `chr5_12300000_12400000.haplotypes.tsv`
- `chr5_12300000_12400000.haplotype_frequency.tsv`
- `chr5_12300000_12400000.haplotype_counts_wide.tsv`
- `chr5_12300000_12400000.sample_haplotypes.tsv`
- `chr5_12300000_12400000.variants.tsv`
- intermediate regional/phased VCF files if `keep_intermediate = TRUE`

## Command Line

After installation:

```bash
Rscript inst/scripts/plot_hapnet.R \
  --vcf=example.vcf.gz \
  --region=chr5:12300000-12400000 \
  --metadata=metadata.txt \
  --group_col=population \
  --outdir=hapnet_chr5_geneA \
  --threads=8
```

## Notes

- The default `snp_only = TRUE` is recommended for haplotype networks.
- Unphased genotypes are not accepted unless `allow_unphased = TRUE`; for
  publication analysis, keep the default Beagle phase step.
- Missing alleles are encoded as `N` in the network input.
- Node size is scaled by haplotype count, and pie slices show group composition.
