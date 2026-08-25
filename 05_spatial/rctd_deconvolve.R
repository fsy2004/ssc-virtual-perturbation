# =============================================================================
# rctd_deconvolve.R  ·  RCTD (spacexr) reference-based deconvolution helper
# -----------------------------------------------------------------------------
# Called as a SUBPROCESS by 06_spatial_immune_coupling.py (Python cannot run
# spacexr). This is a THIN adaptation of reusable module 505_spatial_advanced.R
# (uses the SAME real spacexr API: Reference / SpatialRNA / create.RCTD /
# run.RCTD, doublet_mode='full'). It does NOT invent any API.
#
# Env on server: a dedicated R env with spacexr (module 505 README install
#   recipe: gh-proxy clone dmcable/spacexr + remotes::install_local). Note the
#   env name in the Python config (R_ENV). If spacexr is missing the Python
#   side catches a non-zero exit and falls back to a marker-score proxy.
#
# I/O contract (all paths passed as args; NOTHING hardcoded):
#   Rscript rctd_deconvolve.R <ref_mtx> <ref_genes> <ref_barcodes> <ref_labels> \
#                             <sp_mtx> <sp_genes> <sp_barcodes> <sp_coords> \
#                             <out_weights_csv> [max_cores]
#   ref_*  : single-cell reference (Tabib/Gur fibroblast+immune atlas subset),
#            MatrixMarket counts (gene x cell), genes/barcodes txt, labels csv
#            (barcode,celltype).
#   sp_*   : one spatial slice, MatrixMarket counts (gene x spot), genes/
#            barcodes txt, coords csv (barcode,x,y).
#   out    : spot x celltype proportion csv (row-normalised).
# Idempotent: if <out_weights_csv> exists and is non-empty, exits 0 immediately.
# Reproducible: set.seed(42).
# =============================================================================
suppressWarnings(suppressMessages({
  library(Matrix)
}))
set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 9) stop("need 9 args (see header)")
ref_mtx <- args[1]; ref_genes <- args[2]; ref_bc <- args[3]; ref_lab <- args[4]
sp_mtx  <- args[5]; sp_genes  <- args[6]; sp_bc  <- args[7]; sp_coord <- args[8]
out_csv <- args[9]
max_cores <- if (length(args) >= 10) as.integer(args[10]) else 4L
# A8: doublet_mode is PLATFORM-dependent (spacexr Visium worked-example uses 'doublet').
# Passed as arg[11] by the Python caller: 55um Visium/visium_mtx/Zenodo-Visium -> 'doublet'
# (each ~10um spot is 1-2 cells), Stereo-seq/Xenium -> 'full'. Default 'doublet'
# (Stereo/Xenium should be passed 'full' explicitly).
doublet_mode <- if (length(args) >= 11) args[11] else "doublet"
if (!doublet_mode %in% c("doublet", "full", "multi")) doublet_mode <- "doublet"

# --- idempotent skip ---------------------------------------------------------
if (file.exists(out_csv) && file.info(out_csv)$size > 10) {
  cat(sprintf("[RCTD] %s exists -> skip\n", basename(out_csv))); quit(status = 0)
}
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)

# spacexr must be present; a missing package -> non-zero exit -> Python fallback
ok <- requireNamespace("spacexr", quietly = TRUE)
if (!ok) { cat("[RCTD] spacexr NOT installed -> abort (Python will fall back)\n"); quit(status = 3) }
suppressWarnings(suppressMessages(library(spacexr)))

readmm <- function(mtx, genes, bc) {
  m <- as(Matrix::readMM(mtx), "CsparseMatrix")        # gene x cell/spot
  rownames(m) <- readLines(genes); colnames(m) <- readLines(bc)
  m
}

cat("[RCTD] loading reference ...\n")
ref_counts <- readmm(ref_mtx, ref_genes, ref_bc)
lab <- read.csv(ref_lab, stringsAsFactors = FALSE)      # barcode,celltype
ct  <- factor(setNames(lab$celltype, lab$barcode)[colnames(ref_counts)])
# spacexr requires >= a few cells/type; drop tiny classes to avoid run.RCTD abort
keep_ct <- names(which(table(ct) >= 25))
sel <- ct %in% keep_ct
ref_counts <- ref_counts[, sel]; ct <- droplevels(ct[sel])
cat(sprintf("[RCTD] reference %d genes x %d cells, %d types: %s\n",
            nrow(ref_counts), ncol(ref_counts), nlevels(ct),
            paste(levels(ct), collapse = ", ")))
reference <- Reference(ref_counts, ct)

cat("[RCTD] loading spatial slice ...\n")
sp_counts <- readmm(sp_mtx, sp_genes, sp_bc)
coords <- read.csv(sp_coord, row.names = 1)             # barcode,x,y
coords <- coords[colnames(sp_counts), c("x", "y"), drop = FALSE]
puck <- SpatialRNA(coords, as(sp_counts, "CsparseMatrix"))

cat(sprintf("[RCTD] create.RCTD (max_cores=%d, doublet_mode=%s) ...\n", max_cores, doublet_mode))
myRCTD <- create.RCTD(puck, reference, max_cores = max_cores, CELL_MIN_INSTANCE = 20)
# A8: doublet_mode per platform. 'doublet' = spacexr Visium worked-example (each ~10um
# Visium spot resolves 1-2 cells); 'full' = every spot a full weighted mixture (Stereo/Xenium).
myRCTD <- run.RCTD(myRCTD, doublet_mode = doublet_mode)
# Extract per-spot cell-type PROPORTION weights (spot × n_celltypes). Both full and doublet
# modes populate @results$weights with the N×n_celltypes proportion matrix (gather_results).
# ★2026-07-07 修: doublet 分支原误优先 @results$weights_doublet(N×2 first/second_type 混合权重),
# 破坏 spot×celltype 契约(06 下游只读到 2 列) -> 改优先 @results$weights, weights_doublet 仅兜底。
if (doublet_mode == "full") {
  W <- as.matrix(myRCTD@results$weights)
} else {
  res <- myRCTD@results
  if (!is.null(res$weights)) {
    W <- as.matrix(res$weights)
  } else if (!is.null(res$weights_doublet)) {
    W <- as.matrix(res$weights_doublet)
  } else {
    stop("[RCTD] doublet_mode results have no weights/weights_doublet matrix")
  }
}
W <- sweep(W, 1, rowSums(W) + 1e-9, "/")
write.csv(W, out_csv)
cat(sprintf("[RCTD] wrote %s (%d spots x %d types)\n", basename(out_csv), nrow(W), ncol(W)))
