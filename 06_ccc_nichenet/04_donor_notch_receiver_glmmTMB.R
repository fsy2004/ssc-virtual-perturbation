suppressPackageStartupMessages({
  library(glmmTMB)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("usage: Rscript 04_donor_notch_receiver_glmmTMB.R counts.csv outdir thresholds_csv")
}
counts_path <- args[[1]]
outdir <- args[[2]]
thresholds <- sort(unique(as.integer(strsplit(args[[3]], ",", fixed = TRUE)[[1]])))
if (!length(thresholds) || any(!is.finite(thresholds)) || any(thresholds < 1)) {
  stop("thresholds must be positive integers")
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dat <- read.csv(counts_path, stringsAsFactors = FALSE)
dat <- dat[dat$condition %in% c("HC", "SSc"), ]
dat$condition <- factor(dat$condition, levels = c("HC", "SSc"))
dat$cohort <- factor(dat$cohort)
genes <- c("NOTCH1", "NOTCH2", "NOTCH3", "HES1")

fit_one <- function(d) {
  d <- droplevels(d)
  design <- model.matrix(~ condition + cohort, data = d)
  rank_deficient <- qr(design)$rank < ncol(design)
  result <- list(
    n_donors = length(unique(d$donor_id)),
    n_SSc = length(unique(d$donor_id[d$condition == "SSc"])),
    n_HC = length(unique(d$donor_id[d$condition == "HC"])),
    n_cohorts = nlevels(d$cohort),
    rank_deficient = rank_deficient,
    converged = NA,
    pdHess = NA,
    dispersion = NA,
    estimate_log_or = NA,
    se_log_or = NA,
    OR = NA,
    CI_low = NA,
    CI_high = NA,
    p_value = NA,
    model_message = NA
  )
  if (result$n_SSc < 2 || result$n_HC < 2 || any(d$total <= 0)) {
    result$model_message <- "insufficient nonzero donor denominators"
    return(result)
  }
  if (rank_deficient) {
    result$model_message <- "rank-deficient condition-plus-cohort design"
    return(result)
  }
  fit <- tryCatch(
    glmmTMB(
      cbind(positive, total - positive) ~ condition + cohort,
      family = betabinomial(link = "logit"),
      data = d
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    result$model_message <- conditionMessage(fit)
    return(result)
  }
  coefficient_table <- summary(fit)$coefficients$cond
  if (!"conditionSSc" %in% rownames(coefficient_table)) {
    result$model_message <- "conditionSSc coefficient absent"
    return(result)
  }
  estimate <- coefficient_table["conditionSSc", "Estimate"]
  standard_error <- coefficient_table["conditionSSc", "Std. Error"]
  p_value <- coefficient_table["conditionSSc", "Pr(>|z|)"]
  pearson <- residuals(fit, type = "pearson")
  residual_df <- df.residual(fit)
  result$converged <- isTRUE(fit$sdr$pdHess) && fit$fit$convergence == 0
  result$pdHess <- isTRUE(fit$sdr$pdHess)
  result$dispersion <- ifelse(
    is.finite(residual_df) && residual_df > 0,
    sum(pearson^2) / residual_df,
    NA
  )
  result$estimate_log_or <- unname(estimate)
  result$se_log_or <- unname(standard_error)
  result$OR <- unname(exp(estimate))
  result$CI_low <- unname(exp(estimate - 1.96 * standard_error))
  result$CI_high <- unname(exp(estimate + 1.96 * standard_error))
  result$p_value <- unname(p_value)
  result$model_message <- ifelse(result$converged, "ok", paste("convergence", fit$fit$convergence))
  result
}

fit_frame <- function(d, threshold, omitted_cohort = "none") {
  answer <- fit_one(d)
  data.frame(
    threshold_min_total = threshold,
    gene = d$gene[1],
    omitted_cohort = omitted_cohort,
    n_donors = answer$n_donors,
    n_SSc = answer$n_SSc,
    n_HC = answer$n_HC,
    n_cohorts = answer$n_cohorts,
    OR = answer$OR,
    CI_low = answer$CI_low,
    CI_high = answer$CI_high,
    p_value = answer$p_value,
    estimate_log_or = answer$estimate_log_or,
    se_log_or = answer$se_log_or,
    converged = answer$converged,
    pdHess = answer$pdHess,
    rank_deficient = answer$rank_deficient,
    dispersion_pearson = answer$dispersion,
    model_message = answer$model_message,
    stringsAsFactors = FALSE
  )
}

rows <- list()
diagnostics <- list()
for (threshold in thresholds) {
  for (gene_name in genes) {
    gene_data <- dat[dat$gene == gene_name & dat$total >= threshold, ]
    rows[[length(rows) + 1]] <- fit_frame(gene_data, threshold)
    diagnostics[[length(diagnostics) + 1]] <- data.frame(
      threshold_min_total = threshold,
      gene = gene_name,
      n_donors = nrow(gene_data),
      min_total = ifelse(nrow(gene_data), min(gene_data$total), NA),
      q01_total = ifelse(nrow(gene_data), unname(quantile(gene_data$total, 0.01)), NA),
      median_total = ifelse(nrow(gene_data), median(gene_data$total), NA),
      q99_total = ifelse(nrow(gene_data), unname(quantile(gene_data$total, 0.99)), NA),
      max_total = ifelse(nrow(gene_data), max(gene_data$total), NA),
      zero_positive_donors = ifelse(nrow(gene_data), sum(gene_data$positive == 0), NA),
      all_positive_donors = ifelse(nrow(gene_data), sum(gene_data$positive == gene_data$total), NA),
      stringsAsFactors = FALSE
    )
  }
}
results <- do.call(rbind, rows)
results$q_BH <- NA_real_
results$p_Holm <- NA_real_
for (threshold in thresholds) {
  indices <- which(results$threshold_min_total == threshold & is.finite(results$p_value))
  if (length(indices)) {
    results$q_BH[indices] <- p.adjust(results$p_value[indices], method = "BH")
    results$p_Holm[indices] <- p.adjust(results$p_value[indices], method = "holm")
  }
}
write.csv(results, file.path(outdir, "Figure6d_glmmTMB_results.csv"), row.names = FALSE)
write.csv(do.call(rbind, diagnostics), file.path(outdir, "Figure6d_model_diagnostics.csv"), row.names = FALSE)

primary <- dat[dat$total >= 1, ]
cohort_summary <- aggregate(
  cbind(positive, total) ~ gene + cohort + condition,
  data = primary,
  FUN = sum
)
cohort_summary$fraction <- cohort_summary$positive / cohort_summary$total
write.csv(cohort_summary, file.path(outdir, "Figure6d_cohort_heterogeneity.csv"), row.names = FALSE)

loco_rows <- list()
for (gene_name in genes) {
  gene_data <- primary[primary$gene == gene_name, ]
  for (cohort_name in levels(droplevels(gene_data$cohort))) {
    reduced <- gene_data[gene_data$cohort != cohort_name, ]
    if (nrow(reduced)) {
      loco_rows[[length(loco_rows) + 1]] <- fit_frame(reduced, 1, cohort_name)
    }
  }
}
write.csv(
  do.call(rbind, loco_rows),
  file.path(outdir, "Figure6d_leave_one_cohort_out.csv"),
  row.names = FALSE
)
