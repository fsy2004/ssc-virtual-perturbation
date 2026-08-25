#!/usr/bin/env Rscript
# ============================================================================
# SKIN base-GRN build, STEP 1 (R): merged cell-by-peak matrix + Cicero
# Follows CellOracle official scATAC/Cicero tutorial
#   notebooks/01_ATAC-seq_data_processing/option1_.../01_atacdata_analysis_with_cicero_and_monocle3
# Multi-donor merge via Signac (union peak set + FeatureMatrix re-quantification),
#   as in the Signac "Merging objects" vignette. GSE312129, 14 donors, hg38.
# ============================================================================
suppressMessages({
  library(Signac); library(Seurat); library(GenomicRanges); library(Matrix)
  library(monocle3); library(cicero); library(future)
})
set.seed(2017)
FRAG  <- "/data/ssc/basegrn/frags"
OUT   <- "/data/ssc/basegrn/cicero_out"; dir.create(OUT,  showWarnings=FALSE, recursive=TRUE)
FMDIR <- "/data/ssc/basegrn/fmats";      dir.create(FMDIR, showWarnings=FALSE, recursive=TRUE)
CELLCAP <- 1200L                 # per-donor cap (speed; total ~16k cells; Cicero uses metacells)
options(future.globals.maxSize = 8*1024^3)
plan("multicore", workers = 8)   # FeatureMatrix internal parallelism

# donor label -> c(h5, fragments)
S <- list(
 SSC1 =c("GSM9338143_SSC1_filtered_feature_bc_matrix.h5","GSM9338144_SSC1_atac_fragments.tsv.gz"),
 SSC2 =c("GSM9338145_SSC2_filtered_feature_bc_matrix.h5","GSM9338146_SSC2_atac_fragments.tsv.gz"),
 SSC3 =c("GSM9338147_SSC3_filtered_feature_bc_matrix.h5","GSM9338148_SSC3_atac_fragments.tsv.gz"),
 SSC4 =c("GSM9338149_SSC4_filtered_feature_bc_matrix.h5","GSM9338150_SSC4_atac_fragments.tsv.gz"),
 SSC5 =c("GSM9338151_SSC5_filtered_feature_bc_matrix.h5","GSM9338152_SSC5_atac_fragments.tsv.gz"),
 SSC6 =c("GSM9338153_SSC6_filtered_feature_bc_matrix.h5","GSM9338154_SSC6_atac_fragments.tsv.gz"),
 SSC7 =c("GSM9338155_SSC7_filtered_feature_bc_matrix.h5","GSM9338156_SSC7_atac_fragments.tsv.gz"),
 SSC8 =c("GSM9338157_SSC8_filtered_feature_bc_matrix.h5","GSM9338158_SSC8_atac_fragments.tsv.gz"),
 SSC9 =c("GSM9338159_SSC9_filtered_feature_bc_matrix.h5","GSM9338160_SSC9_atac_fragments.tsv.gz"),
 SSC10=c("GSM9338161_SSC10_filtered_feature_bc_matrix.h5","GSM9338162_SSC10_atac_fragments.tsv.gz"),
 HC1  =c("GSM9338163_HC1_filtered_feature_bc_matrix.h5","GSM9338164_HC1_atac_fragments.tsv.gz"),
 HC2  =c("GSM9338165_HC2_filtered_feature_bc_matrix.h5","GSM9338166_HC2_atac_fragments.tsv.gz"),
 HC3  =c("GSM9338167_HC3_filtered_feature_bc_matrix.h5","GSM9338168_HC3_atac_fragments.tsv.gz"),
 HC4  =c("GSM9338169_HC4_filtered_feature_bc_matrix.h5","GSM9338170_HC4_atac_fragments.tsv.gz"))

std_chroms <- paste0("chr", c(1:22,"X","Y"))

# ---- 1. per-donor peaks + capped valid barcodes (from multiome h5) ----------
cat("[1] reading per-donor peaks + barcodes\n"); flush.console()
peak_grs <- list(); cells_list <- list()
for (s in names(S)) {
  cnt <- Read10X_h5(file.path(FRAG, S[[s]][1]))     # list: Gene Expression + Peaks
  pk  <- cnt[["Peaks"]]
  gr  <- StringToGRanges(rownames(pk), sep = c(":","-"))
  gr  <- gr[as.character(seqnames(gr)) %in% std_chroms]
  peak_grs[[s]] <- gr
  cb <- colnames(pk)
  if (length(cb) > CELLCAP) cb <- sample(cb, CELLCAP)
  cells_list[[s]] <- cb
  cat(sprintf("   %-6s peaks=%d cells=%d\n", s, length(gr), length(cb))); flush.console()
}

# ---- 2. union peak set (reduce) --------------------------------------------
combined <- reduce(do.call(c, unname(peak_grs)))
w <- width(combined); combined <- combined[w >= 50 & w <= 5000]
combined <- combined[as.character(seqnames(combined)) %in% std_chroms]
saveRDS(combined, file.path(OUT,"combined_peaks.rds"))
cat(sprintf("[2] combined peak set: %d peaks (width 50-5000)\n", length(combined))); flush.console()

# ---- 3. per-donor FeatureMatrix re-quantification (checkpointed) -----------
for (s in names(S)) {
  rds <- file.path(FMDIR, paste0(s,".rds"))
  if (file.exists(rds)) { cat("   [skip]", s, "\n"); next }
  fragfile <- file.path(FRAG, S[[s]][2])
  frags <- CreateFragmentObject(path=fragfile, cells=cells_list[[s]], validate.fragments=FALSE)
  fm <- FeatureMatrix(fragments=frags, features=combined, cells=cells_list[[s]], process_n=6000)
  if (sum(fm) == 0) stop(sprintf("FeatureMatrix all-zero for %s -> barcode mismatch; ABORT", s))
  colnames(fm) <- paste0(s, "_", colnames(fm))
  saveRDS(fm, rds)
  cat(sprintf("   [fm] %-6s dim=%dx%d nnz=%d\n", s, nrow(fm), ncol(fm), length(fm@x))); flush.console()
}

# ---- 4. merge, binarize, name peaks chr_start_end --------------------------
fmats <- lapply(names(S), function(s) readRDS(file.path(FMDIR, paste0(s,".rds"))))
mat <- do.call(cbind, fmats)                       # rows identical order = combined
mat@x[mat@x > 0] <- 1
site <- paste(as.character(seqnames(combined)), start(combined), end(combined), sep="_")
rownames(mat) <- site
cat(sprintf("[4] merged matrix: %d peaks x %d cells\n", nrow(mat), ncol(mat))); flush.console()

# ---- 5. monocle3 CDS + light QC (tutorial style) ---------------------------
peakinfo <- data.frame(site_name = rownames(mat)); rownames(peakinfo) <- peakinfo$site_name
cellinfo <- data.frame(cells = colnames(mat));     rownames(cellinfo) <- cellinfo$cells
input_cds <- suppressWarnings(new_cell_data_set(mat, cell_metadata=cellinfo, gene_metadata=peakinfo))
input_cds <- monocle3::detect_genes(input_cds)
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0, ]           # drop empty peaks
ncell <- ncol(input_cds)
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) >= ceiling(0.0025*ncell), ]  # mild prevalence (>=0.25% cells)
input_cds <- input_cds[, Matrix::colSums(exprs(input_cds)) >= 200]         # drop low-coverage cells
cat(sprintf("[5] after QC: %d peaks x %d cells\n", nrow(input_cds), ncol(input_cds))); flush.console()

set.seed(2017)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method="LSI")
input_cds <- reduce_dimension(input_cds, reduction_method="UMAP", preprocess_method="LSI")
umap_coords <- reducedDims(input_cds)$UMAP
cicero_cds  <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

# hg38 UCSC main chromosome lengths (chr1-22,X,Y)
chrlen <- data.frame(
  V1=c(paste0("chr",1:22),"chrX","chrY"),
  V2=c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,
       145138636,138394717,133797422,135086622,133275309,114364328,107043718,
       101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,
       156040895,57227415))
cat("[6] running Cicero ...\n"); flush.console()
conns <- run_cicero(cicero_cds, chrlen)

all_peaks <- row.names(exprs(input_cds))
write.csv(all_peaks, file.path(OUT,"all_peaks.csv"))
write.csv(conns,     file.path(OUT,"cicero_connections.csv"))
cat(sprintf("DONE_R  n_peaks=%d  n_conns=%d\n", length(all_peaks), nrow(conns)))
