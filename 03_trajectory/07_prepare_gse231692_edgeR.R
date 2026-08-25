#!/usr/bin/env Rscript

# Prepare the GSE231692 baseline SSc/HC RNA-seq matrix with the published
# edgeR workflow before CollecTRI/ULM activity inference.

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(edgeR)
  library(org.Hs.eg.db)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: 07_prepare_gse231692_edgeR.R ",
    "<GSE231692_count.tsv.gz> <candidate_sample_metadata.csv> <output_dir>"
  )
}

count_path <- normalizePath(args[[1]], mustWork = TRUE)
metadata_path <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- normalizePath(args[[3]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
metadata <- metadata[metadata$accession == "GSE231692", , drop = FALSE]
is_hc <- metadata$condition_screen == "HC"
is_baseline_ssc <- metadata$condition_screen == "SSc" &
  metadata$timepoint_screen == "baseline"
metadata$analysis_keep <- is_hc | is_baseline_ssc
metadata$exclusion_reason <- ifelse(
  metadata$analysis_keep,
  "",
  "post-baseline SSc biopsy"
)
metadata <- metadata[order(metadata$geo_accession), , drop = FALSE]

if (sum(is_hc) != 14L || sum(is_baseline_ssc) != 36L) {
  stop(
    "Unexpected GSE231692 baseline composition: ",
    sum(is_baseline_ssc), " SSc and ", sum(is_hc), " HC"
  )
}
if (anyDuplicated(metadata$title) || anyDuplicated(metadata$geo_accession)) {
  stop("GSE231692 sample titles or GEO accessions are not unique")
}

counts <- read.delim(
  gzfile(count_path),
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (any(counts < 0, na.rm = TRUE) || any(abs(counts - round(counts)) > 1e-8, na.rm = TRUE)) {
  stop("GSE231692 matrix does not contain non-negative integer counts")
}

kept <- metadata[metadata$analysis_keep, , drop = FALSE]
missing_titles <- setdiff(kept$title, colnames(counts))
if (length(missing_titles)) {
  stop("Retained samples absent from count matrix: ", paste(missing_titles, collapse = ", "))
}
counts <- as.matrix(counts[, kept$title, drop = FALSE])
storage.mode(counts) <- "integer"
colnames(counts) <- kept$geo_accession

ensembl <- sub("\\..*$", "", rownames(counts))
symbols <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = unique(ensembl),
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)
symbol_by_row <- unname(symbols[ensembl])
mapped <- !is.na(symbol_by_row) & nzchar(symbol_by_row)
counts_by_symbol <- rowsum(
  counts[mapped, , drop = FALSE],
  group = symbol_by_row[mapped],
  reorder = FALSE
)

group <- factor(kept$condition_screen, levels = c("HC", "SSc"))
y <- DGEList(counts = counts_by_symbol, group = group)
expressed <- filterByExpr(y, group = group)
y <- y[expressed, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")
log_cpm <- cpm(y, log = TRUE, prior.count = 2)

expression_out <- file.path(output_dir, "GSE231692_edgeR_logCPM.csv")
sample_qc_out <- file.path(output_dir, "GSE231692_edgeR_sample_qc.csv")
gene_qc_out <- file.path(output_dir, "GSE231692_edgeR_gene_qc.csv")
mapping_out <- file.path(output_dir, "GSE231692_edgeR_sample_selection.csv")

write.csv(t(log_cpm), expression_out, quote = FALSE)
write.csv(
  data.frame(
    geo_accession = colnames(y),
    subject_id = kept$subject_id[match(colnames(y), kept$geo_accession)],
    condition = kept$condition_screen[match(colnames(y), kept$geo_accession)],
    library_size = y$samples$lib.size,
    TMM_norm_factor = y$samples$norm.factors,
    effective_library_size = y$samples$lib.size * y$samples$norm.factors,
    check.names = FALSE
  ),
  sample_qc_out,
  row.names = FALSE,
  quote = FALSE
)
write.csv(
  data.frame(
    metric = c(
      "source_ensembl_rows",
      "mapped_ensembl_rows",
      "unique_symbols_before_filterByExpr",
      "symbols_after_filterByExpr",
      "retained_samples",
      "retained_SSc",
      "retained_HC"
    ),
    value = c(
      nrow(counts),
      sum(mapped),
      nrow(counts_by_symbol),
      nrow(y),
      ncol(y),
      sum(group == "SSc"),
      sum(group == "HC")
    )
  ),
  gene_qc_out,
  row.names = FALSE,
  quote = FALSE
)
write.csv(
  metadata[, c(
    "geo_accession",
    "title",
    "subject_id",
    "condition_screen",
    "timepoint_screen",
    "analysis_keep",
    "exclusion_reason"
  )],
  mapping_out,
  row.names = FALSE,
  quote = FALSE
)

cat(
  sprintf(
    "GSE231692: %d baseline subjects (%d SSc, %d HC), %d symbols after filterByExpr\n",
    ncol(y), sum(group == "SSc"), sum(group == "HC"), nrow(y)
  )
)
