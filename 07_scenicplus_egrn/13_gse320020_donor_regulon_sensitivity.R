#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(edgeR)
  library(limma)
  library(decoupleR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "usage: Rscript 13_gse320020_donor_regulon_sensitivity.R ",
    "<RData.gz> <approved_source_mapping.tsv> <collectri.csv> <output_dir>"
  )
}

input_path <- normalizePath(args[[1]], mustWork = TRUE)
mapping_path <- normalizePath(args[[2]], mustWork = TRUE)
network_path <- normalizePath(args[[3]], mustWork = TRUE)
output_dir <- args[[4]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

focal_tfs <- c("HES1", "SMAD3", "FOSB", "MEF2C")
focal_genes <- c(focal_tfs, "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4")
condition_colors <- c(HC = "#0072B2", SSc = "#D55E00")

as_flag <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes")
}

write_tsv <- function(x, path) {
  write.table(
    x,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}

write_tsv_gz <- function(x, path) {
  connection <- gzfile(path, open = "wt")
  on.exit(close(connection), add = TRUE)
  write_tsv(x, connection)
}

mapping <- read.delim(
  mapping_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
required_mapping <- c(
  "object_name",
  "raw_donor",
  "raw_condition",
  "raw_source",
  "n_cells",
  "approved_donor_id",
  "approved_condition",
  "source_study",
  "include",
  "approved_by"
)
missing_mapping <- setdiff(required_mapping, colnames(mapping))
if (length(missing_mapping)) {
  stop("Approved mapping lacks columns: ", paste(missing_mapping, collapse = ", "))
}
mapping$include <- as_flag(mapping$include)
included_mapping <- mapping[mapping$include, , drop = FALSE]
if (!nrow(included_mapping)) {
  stop("Approved mapping retains no rows")
}
if (anyDuplicated(included_mapping$approved_donor_id)) {
  stop("Included donor identifiers must be unique")
}
if (!all(included_mapping$approved_condition %in% c("HC", "SSc"))) {
  stop("Included conditions must be HC or SSc")
}
if (!all(included_mapping$object_name == "ASCT9V9_Fibro")) {
  stop("Only the approved ASCT9V9_Fibro object is allowed")
}

loaded_env <- new.env(parent = emptyenv())
loaded_names <- load(gzfile(input_path), envir = loaded_env)
if (!"ASCT9V9_Fibro" %in% loaded_names) {
  stop("ASCT9V9_Fibro is absent from the RData archive")
}
object <- get("ASCT9V9_Fibro", envir = loaded_env)
if (!inherits(object, "Seurat")) {
  stop("ASCT9V9_Fibro is not a Seurat object")
}

metadata <- object[[]]
metadata$cell_id <- rownames(metadata)
required_metadata <- c("library_id", "health", "seq")
missing_metadata <- setdiff(required_metadata, colnames(metadata))
if (length(missing_metadata)) {
  stop("Fibroblast metadata lacks columns: ", paste(missing_metadata, collapse = ", "))
}

mapping$key <- paste(
  mapping$object_name,
  mapping$raw_donor,
  mapping$raw_condition,
  mapping$raw_source,
  sep = "\r"
)
metadata$key <- paste(
  "ASCT9V9_Fibro",
  as.character(metadata$library_id),
  as.character(metadata$health),
  as.character(metadata$seq),
  sep = "\r"
)
map_index <- match(metadata$key, mapping$key)
if (anyNA(map_index)) {
  stop("At least one fibroblast cell has no approved donor/source mapping")
}
cell_mapping <- mapping[map_index, , drop = FALSE]
keep_cells <- cell_mapping$include
if (!any(keep_cells)) {
  stop("No fibroblast cells remain after approved mapping")
}

cell_metadata <- data.frame(
  cell_id = metadata$cell_id[keep_cells],
  donor_id = cell_mapping$approved_donor_id[keep_cells],
  condition = cell_mapping$approved_condition[keep_cells],
  source_study = cell_mapping$source_study[keep_cells],
  raw_assay = cell_mapping$raw_source[keep_cells],
  stringsAsFactors = FALSE
)
observed_counts <- as.data.frame(table(cell_metadata$donor_id))
colnames(observed_counts) <- c("approved_donor_id", "observed_n_cells")
observed_counts$approved_donor_id <- as.character(
  observed_counts$approved_donor_id
)
donor_metadata <- merge(
  included_mapping,
  observed_counts,
  by = "approved_donor_id",
  all.x = TRUE,
  sort = FALSE
)
if (any(donor_metadata$n_cells != donor_metadata$observed_n_cells)) {
  bad <- donor_metadata[
    donor_metadata$n_cells != donor_metadata$observed_n_cells,
    c("approved_donor_id", "n_cells", "observed_n_cells"),
    drop = FALSE
  ]
  stop(
    "Approved cell counts differ from the live Seurat object: ",
    paste(capture.output(print(bad, row.names = FALSE)), collapse = " ")
  )
}

donor_metadata <- donor_metadata[
  order(
    match(donor_metadata$approved_condition, c("HC", "SSc")),
    donor_metadata$source_study,
    donor_metadata$approved_donor_id
  ),
  ,
  drop = FALSE
]
donor_ids <- donor_metadata$approved_donor_id
donor_metadata$analysis_same_study <-
  donor_metadata$source_study == "GSE320020"
donor_metadata$analysis_integrated <- TRUE
write_tsv(
  donor_metadata[
    ,
    c(
      "approved_donor_id",
      "approved_condition",
      "source_study",
      "raw_source",
      "observed_n_cells",
      "analysis_same_study",
      "analysis_integrated",
      "approval_note"
    )
  ],
  file.path(output_dir, "GSE320020_donor_metadata.tsv")
)

counts <- LayerData(object[["RNA"]], layer = "counts")
if (!inherits(counts, "sparseMatrix")) {
  counts <- as(counts, "dgCMatrix")
}
cell_order <- match(cell_metadata$cell_id, colnames(counts))
if (anyNA(cell_order)) {
  stop("Approved fibroblast cells are absent from the RNA counts layer")
}
counts <- counts[, cell_order, drop = FALSE]
membership <- sparseMatrix(
  i = seq_len(nrow(cell_metadata)),
  j = match(cell_metadata$donor_id, donor_ids),
  x = 1,
  dims = c(nrow(cell_metadata), length(donor_ids)),
  dimnames = list(cell_metadata$cell_id, donor_ids)
)
pseudobulk <- counts %*% membership
storage.mode(pseudobulk@x) <- "double"
if (!identical(colnames(pseudobulk), donor_ids)) {
  stop("Pseudobulk donor order is inconsistent with donor metadata")
}
if (any(Matrix::colSums(pseudobulk) <= 0)) {
  stop("At least one donor has a zero pseudobulk library")
}

pseudobulk_frame <- data.frame(
  gene = rownames(pseudobulk),
  as.matrix(pseudobulk),
  check.names = FALSE
)
write_tsv_gz(
  pseudobulk_frame,
  file.path(output_dir, "GSE320020_donor_pseudobulk_counts.tsv.gz")
)

dge_all <- DGEList(counts = pseudobulk)
dge_all <- calcNormFactors(dge_all, method = "TMM")
log_cpm <- cpm(dge_all, log = TRUE, prior.count = 2)
log_cpm_frame <- data.frame(
  gene = rownames(log_cpm),
  log_cpm,
  check.names = FALSE
)
write_tsv_gz(
  log_cpm_frame,
  file.path(output_dir, "GSE320020_donor_TMM_logCPM.tsv.gz")
)

network <- read.csv(network_path, stringsAsFactors = FALSE)
required_network <- c("source", "target", "weight")
missing_network <- setdiff(required_network, colnames(network))
if (length(missing_network)) {
  stop("CollecTRI network lacks columns: ", paste(missing_network, collapse = ", "))
}
network <- network[
  is.finite(network$weight) &
    network$target %in% rownames(log_cpm),
  required_network,
  drop = FALSE
]
if (!"HES1" %in% network$source) {
  stop("HES1 is absent after intersecting CollecTRI with expressed genes")
}

ulm_long <- run_ulm(
  mat = log_cpm,
  network = network,
  .source = source,
  .target = target,
  .mor = weight,
  sparse = FALSE,
  center = FALSE,
  minsize = 5L
)
score_column <- if ("score" %in% colnames(ulm_long)) {
  "score"
} else if ("estimate" %in% colnames(ulm_long)) {
  "estimate"
} else {
  NA_character_
}
sample_column <- if ("condition" %in% colnames(ulm_long)) {
  "condition"
} else if ("sample" %in% colnames(ulm_long)) {
  "sample"
} else if ("sample_id" %in% colnames(ulm_long)) {
  "sample_id"
} else {
  NA_character_
}
if (is.na(score_column) || is.na(sample_column)) {
  stop(
    "Unexpected decoupleR ULM columns: ",
    paste(colnames(ulm_long), collapse = ", ")
  )
}
if (!"source" %in% colnames(ulm_long)) {
  stop("decoupleR ULM output lacks the source column")
}
ulm_export <- data.frame(
  donor_id = as.character(ulm_long[[sample_column]]),
  TF = as.character(ulm_long$source),
  score = as.numeric(ulm_long[[score_column]]),
  stringsAsFactors = FALSE
)
if (any(!ulm_export$donor_id %in% donor_ids)) {
  stop("decoupleR returned unknown donor identifiers")
}
write_tsv(
  ulm_export,
  file.path(output_dir, "GSE320020_collectri_ULM_donor_activity_long.tsv")
)

tf_levels <- sort(unique(ulm_export$TF))
activity <- matrix(
  NA_real_,
  nrow = length(tf_levels),
  ncol = length(donor_ids),
  dimnames = list(tf_levels, donor_ids)
)
activity[
  cbind(
    match(ulm_export$TF, tf_levels),
    match(ulm_export$donor_id, donor_ids)
  )
] <- ulm_export$score
if (anyNA(activity)) {
  stop("The donor-by-TF ULM matrix is incomplete")
}

activity_export <- data.frame(
  TF = rownames(activity),
  activity,
  check.names = FALSE
)
write_tsv(
  activity_export,
  file.path(output_dir, "GSE320020_collectri_ULM_donor_activity.tsv")
)

limma_table <- function(matrix_values, selected_donors, feature_name) {
  selected_metadata <- donor_metadata[
    match(selected_donors, donor_metadata$approved_donor_id),
    ,
    drop = FALSE
  ]
  group <- factor(
    selected_metadata$approved_condition,
    levels = c("HC", "SSc")
  )
  if (sum(group == "HC") < 2L || sum(group == "SSc") < 2L) {
    stop(feature_name, " design has fewer than two donors in a group")
  }
  design <- model.matrix(~ group)
  fit <- eBayes(
    lmFit(matrix_values[, selected_donors, drop = FALSE], design),
    robust = TRUE
  )
  result <- topTable(
    fit,
    coef = "groupSSc",
    number = Inf,
    sort.by = "none",
    confint = TRUE
  )
  result$feature <- rownames(result)
  result$mean_HC <- rowMeans(
    matrix_values[, selected_donors[group == "HC"], drop = FALSE]
  )[result$feature]
  result$mean_SSc <- rowMeans(
    matrix_values[, selected_donors[group == "SSc"], drop = FALSE]
  )[result$feature]
  result$n_HC <- sum(group == "HC")
  result$n_SSc <- sum(group == "SSc")
  result$rank_positive_effect <- rank(-result$logFC, ties.method = "min")
  result$rank_absolute_t <- rank(-abs(result$t), ties.method = "min")
  result[
    ,
    c(
      "feature",
      "mean_HC",
      "mean_SSc",
      "logFC",
      "CI.L",
      "CI.R",
      "t",
      "P.Value",
      "adj.P.Val",
      "n_HC",
      "n_SSc",
      "rank_positive_effect",
      "rank_absolute_t"
    )
  ]
}

voom_table <- function(selected_donors, design_label) {
  selected_metadata <- donor_metadata[
    match(selected_donors, donor_metadata$approved_donor_id),
    ,
    drop = FALSE
  ]
  group <- factor(
    selected_metadata$approved_condition,
    levels = c("HC", "SSc")
  )
  design <- model.matrix(~ group)
  dge <- DGEList(counts = pseudobulk[, selected_donors, drop = FALSE])
  keep <- filterByExpr(dge, design = design)
  keep[rownames(dge) %in% focal_genes] <- TRUE
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  voom_fit <- voomWithQualityWeights(dge, design, plot = FALSE)
  fit <- eBayes(
    lmFit(voom_fit, design),
    robust = TRUE,
    trend = FALSE
  )
  result <- topTable(
    fit,
    coef = "groupSSc",
    number = Inf,
    sort.by = "none",
    confint = TRUE
  )
  result$feature <- rownames(result)
  result$mean_HC <- rowMeans(
    voom_fit$E[, group == "HC", drop = FALSE]
  )[result$feature]
  result$mean_SSc <- rowMeans(
    voom_fit$E[, group == "SSc", drop = FALSE]
  )[result$feature]
  result$n_HC <- sum(group == "HC")
  result$n_SSc <- sum(group == "SSc")
  result$rank_positive_effect <- rank(-result$logFC, ties.method = "min")
  result$rank_absolute_t <- rank(-abs(result$t), ties.method = "min")
  sample_weights <- data.frame(
    design = design_label,
    donor_id = selected_donors,
    sample_weight = as.numeric(voom_fit$targets$sample.weights),
    stringsAsFactors = FALSE
  )
  list(
    statistics = result[
      ,
      c(
        "feature",
        "mean_HC",
        "mean_SSc",
        "logFC",
        "CI.L",
        "CI.R",
        "t",
        "P.Value",
        "adj.P.Val",
        "n_HC",
        "n_SSc",
        "rank_positive_effect",
        "rank_absolute_t"
      )
    ],
    sample_weights = sample_weights
  )
}

designs <- list(
  same_study = donor_metadata$approved_donor_id[
    donor_metadata$analysis_same_study
  ],
  integrated = donor_metadata$approved_donor_id
)
regulon_results <- list()
expression_results <- list()
sample_weight_rows <- list()

for (design_name in names(designs)) {
  selected_donors <- designs[[design_name]]
  regulon_results[[design_name]] <- limma_table(
    activity,
    selected_donors,
    paste0(design_name, " regulon")
  )
  expression_fit <- voom_table(selected_donors, design_name)
  expression_results[[design_name]] <- expression_fit$statistics
  sample_weight_rows[[design_name]] <- expression_fit$sample_weights
  write_tsv(
    regulon_results[[design_name]],
    file.path(
      output_dir,
      paste0("GSE320020_regulon_", design_name, "_limma.tsv")
    )
  )
  write_tsv(
    expression_results[[design_name]],
    file.path(
      output_dir,
      paste0("GSE320020_expression_", design_name, "_voom_limma.tsv")
    )
  )
}
write_tsv(
  do.call(rbind, sample_weight_rows),
  file.path(output_dir, "GSE320020_expression_voom_sample_weights.tsv")
)

focal_rows <- list()
for (design_name in names(designs)) {
  regulon <- regulon_results[[design_name]]
  regulon <- regulon[regulon$feature %in% focal_tfs, , drop = FALSE]
  regulon$design <- design_name
  regulon$feature_type <- "CollecTRI_ULM"
  expression <- expression_results[[design_name]]
  expression <- expression[expression$feature %in% focal_genes, , drop = FALSE]
  expression$design <- design_name
  expression$feature_type <- "pseudobulk_expression"
  focal_rows[[paste0(design_name, "_regulon")]] <- regulon
  focal_rows[[paste0(design_name, "_expression")]] <- expression
}
focal <- do.call(rbind, focal_rows)
rownames(focal) <- NULL
focal <- focal[
  ,
  c(
    "design",
    "feature_type",
    "feature",
    "mean_HC",
    "mean_SSc",
    "logFC",
    "CI.L",
    "CI.R",
    "t",
    "P.Value",
    "adj.P.Val",
    "n_HC",
    "n_SSc",
    "rank_positive_effect",
    "rank_absolute_t"
  )
]
write_tsv(
  focal,
  file.path(output_dir, "GSE320020_focal_donor_sensitivity.tsv")
)

draw_candidate_figure <- function(path, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") {
    pdf(path, width = 11, height = 8.5, useDingbats = FALSE)
  } else {
    png(path, width = 3300, height = 2550, res = 300, type = "cairo")
  }
  on.exit(dev.off(), add = TRUE)
  layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE))
  par(
    mar = c(5.0, 5.0, 3.2, 1.2),
    mgp = c(2.8, 0.8, 0),
    tcl = -0.25,
    family = "sans"
  )

  bar_colors <- condition_colors[donor_metadata$approved_condition]
  barplot(
    donor_metadata$observed_n_cells,
    names.arg = donor_metadata$approved_donor_id,
    las = 2,
    col = bar_colors,
    border = NA,
    ylab = "Fibroblasts per donor",
    main = "a  Approved donor and source audit"
  )
  legend(
    "topright",
    legend = names(condition_colors),
    fill = condition_colors,
    border = NA,
    bty = "n",
    horiz = TRUE
  )
  mtext(
    paste0(
      "9 SSc; 2 same-study HC; 7 provenance-mapped external HC"
    ),
    side = 3,
    line = 0.4,
    cex = 0.75
  )

  hes1 <- activity["HES1", donor_metadata$approved_donor_id]
  group_numeric <- ifelse(
    donor_metadata$approved_condition == "HC",
    1,
    2
  )
  plot(
    jitter(group_numeric, amount = 0.08),
    hes1,
    pch = ifelse(donor_metadata$source_study == "GSE320020", 21, 24),
    bg = bar_colors,
    col = "#222222",
    xaxt = "n",
    xlim = c(0.65, 2.35),
    xlab = "",
    ylab = "HES1 CollecTRI/ULM score",
    main = "b  Donor-level HES1 activity"
  )
  axis(1, at = 1:2, labels = c("HC", "SSc"))
  segments(
    0.82,
    mean(hes1[group_numeric == 1]),
    1.18,
    mean(hes1[group_numeric == 1]),
    lwd = 2.5
  )
  segments(
    1.82,
    mean(hes1[group_numeric == 2]),
    2.18,
    mean(hes1[group_numeric == 2]),
    lwd = 2.5
  )
  legend(
    "topright",
    legend = c("GSE320020 donor", "External HC donor"),
    pch = c(21, 24),
    pt.bg = "#999999",
    bty = "n",
    cex = 0.75
  )

  hes1_forest <- focal[
    focal$feature == "HES1" &
      focal$feature_type == "CollecTRI_ULM",
    ,
    drop = FALSE
  ]
  hes1_forest <- hes1_forest[
    match(c("same_study", "integrated"), hes1_forest$design),
    ,
    drop = FALSE
  ]
  y <- rev(seq_len(nrow(hes1_forest)))
  x_range <- range(
    c(hes1_forest$CI.L, hes1_forest$CI.R, 0),
    finite = TRUE
  )
  plot(
    NA,
    xlim = x_range,
    ylim = c(0.5, nrow(hes1_forest) + 0.5),
    yaxt = "n",
    xlab = "Mean SSc-HC activity difference",
    ylab = "",
    main = "c  HES1 donor-level estimates"
  )
  abline(v = 0, col = "#777777", lty = 2)
  segments(
    hes1_forest$CI.L,
    y,
    hes1_forest$CI.R,
    y,
    lwd = 2,
    col = "#444444"
  )
  points(
    hes1_forest$logFC,
    y,
    pch = 21,
    bg = condition_colors[["SSc"]],
    cex = 1.2
  )
  axis(
    2,
    at = y,
    labels = c("Same study: 9 vs 2", "Integrated: 9 vs 9"),
    las = 1
  )

  same <- regulon_results$same_study
  integrated <- regulon_results$integrated
  joined <- merge(
    same[, c("feature", "logFC")],
    integrated[, c("feature", "logFC")],
    by = "feature",
    suffixes = c("_same", "_integrated")
  )
  plot(
    joined$logFC_same,
    joined$logFC_integrated,
    pch = 16,
    cex = 0.45,
    col = "#A0A0A080",
    xlab = "Same-study SSc-HC activity difference",
    ylab = "Integrated SSc-HC activity difference",
    main = "d  Full-regulon effect agreement"
  )
  abline(h = 0, v = 0, col = "#BBBBBB", lty = 2)
  abline(0, 1, col = "#555555", lty = 3)
  highlighted <- joined$feature %in% focal_tfs
  points(
    joined$logFC_same[highlighted],
    joined$logFC_integrated[highlighted],
    pch = 21,
    bg = condition_colors[["SSc"]],
    cex = 1.0
  )
  text(
    joined$logFC_same[highlighted],
    joined$logFC_integrated[highlighted],
    labels = joined$feature[highlighted],
    pos = 4,
    cex = 0.7,
    xpd = NA
  )
}

draw_candidate_figure(
  file.path(output_dir, "GSE320020_donor_sensitivity_candidate.pdf"),
  "pdf"
)
draw_candidate_figure(
  file.path(output_dir, "GSE320020_donor_sensitivity_candidate.png"),
  "png"
)

manifest <- c(
  "dataset=GSE320020",
  "object=ASCT9V9_Fibro",
  paste0("input=", input_path),
  paste0("approved_mapping=", mapping_path),
  paste0("collectri_network=", network_path),
  paste0("n_donors=", nrow(donor_metadata)),
  paste0(
    "same_study=",
    sum(
      donor_metadata$analysis_same_study &
        donor_metadata$approved_condition == "SSc"
    ),
    "_SSc_vs_",
    sum(
      donor_metadata$analysis_same_study &
        donor_metadata$approved_condition == "HC"
    ),
    "_HC"
  ),
  paste0(
    "integrated=",
    sum(donor_metadata$approved_condition == "SSc"),
    "_SSc_vs_",
    sum(donor_metadata$approved_condition == "HC"),
    "_HC"
  ),
  "aggregation=raw_counts_summed_within_approved_donor",
  "normalization=edgeR_TMM_logCPM_prior_count_2",
  "regulon=CollecTRI_ULM_minsize_5",
  "regulon_inference=donor_equal_weight_limma_empirical_Bayes",
  paste0(
    "expression_inference=filterByExpr_plus_prespecified_focal_genes;",
    "voomWithQualityWeights;limma_empirical_Bayes"
  ),
  "same_study_and_integrated_designs_reported_separately",
  "integrated_external_HC_sources=GSE264508;GSE138669;GSE288490",
  "evidence_boundary=domain_shift_sensitivity_not_fully_independent_validation"
)
writeLines(
  manifest,
  file.path(output_dir, "GSE320020_donor_sensitivity_manifest.txt")
)
capture.output(
  sessionInfo(),
  file = file.path(output_dir, "GSE320020_donor_sensitivity_sessionInfo.txt")
)
writeLines(
  format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  file.path(output_dir, "GSE320020_DONOR_SENSITIVITY_DONE")
)

cat(
  sprintf(
    "[ok] GSE320020 donor sensitivity completed: %d donors, %d ULM TFs\n",
    nrow(donor_metadata),
    nrow(activity)
  )
)
