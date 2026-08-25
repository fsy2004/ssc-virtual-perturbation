#!/usr/bin/env Rscript

# Donor-pseudobulk differential accessibility for the condition-aware SSc skin
# SCENIC+ extension. This is a new output arm and does not replace the completed
# shared-source SCENIC+ archive.

suppressPackageStartupMessages({
  library(edgeR)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(
    "usage: Rscript 09_ssc_fibroblast_dar_edgeR.R ",
    "<fibroblast_donor_pseudobulk_dir> <output_dir>"
  )
}
input_dir <- args[[1]]
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

counts_path <- file.path(input_dir, "fibroblast_peak_by_donor.mtx.gz")
regions_path <- file.path(input_dir, "regions.tsv.gz")
samples_path <- file.path(input_dir, "samples.tsv")
for (path in c(counts_path, regions_path, samples_path)) {
  if (!file.exists(path)) stop("missing input: ", path)
}

counts <- readMM(gzfile(counts_path))
regions <- read.delim(
  gzfile(regions_path), header = FALSE, stringsAsFactors = FALSE
)[[1]]
samples <- read.delim(samples_path, stringsAsFactors = FALSE)
if (nrow(counts) != length(regions)) stop("region count does not match matrix rows")
if (ncol(counts) != nrow(samples)) stop("sample count does not match matrix columns")
if (!all(samples$condition %in% c("HC", "SSc"))) stop("unexpected condition label")
if (length(unique(samples$sample_id)) != nrow(samples)) stop("duplicate sample_id")

samples$condition <- factor(samples$condition, levels = c("HC", "SSc"))
counts <- round(as.matrix(counts))
storage.mode(counts) <- "integer"
rownames(counts) <- regions
colnames(counts) <- samples$sample_id

y <- DGEList(counts = counts, samples = samples, group = samples$condition)
keep <- filterByExpr(y, group = samples$condition)
if (sum(keep) < 100) stop("fewer than 100 peaks passed donor-level filterByExpr")
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~ condition, data = samples)
if (qr(design)$rank < ncol(design)) stop("rank-deficient condition design")

y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
test <- glmQLFTest(fit, coef = "conditionSSc")
result <- topTags(test, n = Inf, sort.by = "none")$table
result$region <- rownames(result)
result <- result[, c("region", setdiff(colnames(result), "region"))]
result <- result[order(result$FDR, -abs(result$logFC)), ]

write.csv(
  result,
  file.path(output_dir, "fibroblast_SSc_vs_HC_donor_pseudobulk_DAR.csv"),
  row.names = FALSE
)
write.csv(
  samples,
  file.path(output_dir, "fibroblast_SSc_vs_HC_samples_used.csv"),
  row.names = FALSE
)
writeLines(
  c(
    "unit=donor",
    "contrast=SSc_vs_HC",
    "model=edgeR_quasi_likelihood",
    "design=~condition",
    paste0("n_HC=", sum(samples$condition == "HC")),
    paste0("n_SSc=", sum(samples$condition == "SSc")),
    paste0("peaks_tested=", nrow(result)),
    "selection_rule=report_all_peaks_with_BH_FDR",
    "parameter_selection=outcome_blind"
  ),
  file.path(output_dir, "fibroblast_SSc_vs_HC_DAR_manifest.txt")
)
capture.output(
  sessionInfo(),
  file = file.path(output_dir, "fibroblast_SSc_vs_HC_sessionInfo.txt")
)
