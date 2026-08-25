# =============================================================================
# banksy_domains.R  ·  BANKSY spatial-domain segmentation helper
# -----------------------------------------------------------------------------
# Called as a SUBPROCESS by 06_spatial_immune_coupling.py. THIN adaptation of
# reusable module 541_banksy_spatial_domains.R — SAME real BANKSY API
# (computeBanksy / runBanksyPCA / clusterBanksy, lambda=0 baseline plus a
# prespecified spatial value).
# Stripped of the framework/plotting deps so it is self-contained on the server;
# it ONLY writes spot_domain_assignments.csv (Python reads the last column).
# Does NOT invent API. If Banksy is missing -> non-zero exit -> Python falls back
# to an expression-Leiden domain (== the lambda=0 non-spatial baseline).
#
#   Rscript banksy_domains.R --input slice_long.csv --outdir <dir> \
#     --lambda 0.25 --k_geom 15 --npcs 20
#     input  : long-format CSV: spot,x,y,gene1,...,geneN  (no 'domain' truth col)
#     output : <outdir>/spot_domain_assignments.csv (spot,x,y,baseline,spatial)
# Idempotent: skips if the output already exists. Reproducible: set.seed(42).
# =============================================================================
suppressWarnings(suppressMessages({ library(Matrix) }))
set.seed(42)

# --- tiny arg parser (module 541 uses bio_args; here inline, no framework dep) ---
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- which(args == flag); if (length(i) && i < length(args)) args[i + 1] else default
}
input  <- getarg("--input")
outdir <- getarg("--outdir", ".")
# SSc-skin high-resolution precedent: mean-neighbour k=15, gradient-neighbour
# k=30, lambda=0.25 and 20 PCs (Li et al., Ann Rheum Dis 2025,
# doi:10.1016/j.ard.2025.06.002). The Python caller still selects a
# platform-specific lambda and records all values in the output path.
lambda <- as.numeric(getarg("--lambda", "0.25"))
k_geom <- as.numeric(getarg("--k_geom", "15"))
res    <- as.numeric(getarg("--res", "0.8"))
npcs   <- as.numeric(getarg("--npcs", "20"))
stopifnot(!is.null(input))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

out_csv <- file.path(outdir, "spot_domain_assignments.csv")
if (file.exists(out_csv) && file.info(out_csv)$size > 10) {
  cat(sprintf("[BANKSY] %s exists -> skip\n", out_csv)); quit(status = 0)
}

ok <- requireNamespace("Banksy", quietly = TRUE) &&
      requireNamespace("SpatialExperiment", quietly = TRUE)
if (!ok) { cat("[BANKSY] Banksy/SpatialExperiment NOT installed -> abort (Python fallback)\n"); quit(status = 3) }
suppressWarnings(suppressMessages({
  library(Banksy); library(SpatialExperiment); library(SummarizedExperiment); library(S4Vectors)
}))

# --- read long-format slice -> SpatialExperiment (module 541 assembly) ---------
raw <- read.csv(input, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c("x", "y") %in% colnames(raw)))
spot_id <- if ("spot" %in% colnames(raw)) raw$spot else sprintf("spot%05d", seq_len(nrow(raw)))
meta_cols <- intersect(c("spot", "x", "y", "domain"), colnames(raw))
gene_cols <- setdiff(colnames(raw), meta_cols)
expr_mat  <- t(as.matrix(raw[, gene_cols, drop = FALSE]))       # gene x spot
colnames(expr_mat) <- spot_id; mode(expr_mat) <- "numeric"
coords <- as.matrix(raw[, c("x", "y")]); rownames(coords) <- spot_id
se <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = expr_mat), spatialCoords = coords,
  colData = S4Vectors::DataFrame(row.names = spot_id))
cat(sprintf("[BANKSY] %d genes x %d spots\n", nrow(se), ncol(se)))

# --- normalize BEFORE BANKSY (Bioconductor Banksy vignette requirement) --------
# computeBanksy expects a normalized assay, not raw integer counts.
se <- scuttle::computeLibraryFactors(se)
assay(se, "normcounts") <- scuttle::normalizeCounts(se, log = FALSE)

# --- BANKSY: neighbour-augmented features; lambda=0 plus spatial setting -------
se <- computeBanksy(se, assay_name = "normcounts", compute_agf = TRUE,
                    k_geom = c(k_geom, k_geom * 2), verbose = FALSE)
lams <- sort(unique(c(0, lambda)))                                # 0 = non-spatial baseline
se <- runBanksyPCA(se, use_agf = TRUE, lambda = lams, npcs = npcs, seed = 42)
se <- clusterBanksy(se, use_agf = TRUE, lambda = lams, resolution = res,
                    algo = "leiden", seed = 42)
cn <- clusterNames(se)
base_name <- grep("lam0_", cn, value = TRUE)[1]
if (is.na(base_name)) base_name <- grep("lam0\\b", cn, value = TRUE)[1]
spat_name <- setdiff(cn, base_name)[1]
if (is.na(spat_name))                                    # ★2026-07-07: 仅 lambda=0 -> 无空间列, fail-loud
  stop("[BANKSY] need lambda > 0 for a spatial-domain column")

out <- data.frame(spot = spot_id, x = coords[, 1], y = coords[, 2],
                  baseline = as.character(colData(se)[[base_name]]),
                  spatial  = as.character(colData(se)[[spat_name]]),
                  check.names = FALSE)
write.csv(out, out_csv, row.names = FALSE)
cat(sprintf("[BANKSY] wrote %s (baseline %d dom, spatial %d dom)\n", basename(out_csv),
            length(unique(out$baseline)), length(unique(out$spatial))))
