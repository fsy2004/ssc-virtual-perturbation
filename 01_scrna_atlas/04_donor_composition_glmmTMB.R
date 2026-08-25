suppressPackageStartupMessages({
  library(glmmTMB)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: Rscript 04_donor_composition_glmmTMB.R counts.csv outdir")
counts_path <- args[[1]]
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dat <- read.csv(counts_path, stringsAsFactors = FALSE)
dat <- dat[dat$condition %in% c("HC", "SSc") & dat$total > 0, ]
dat$condition <- factor(dat$condition, levels = c("HC", "SSc"))
dat$cohort <- factor(dat$cohort)

fit_one <- function(d) {
  d <- droplevels(d)
  design <- model.matrix(~ condition + cohort, data = d)
  rank_deficient <- qr(design)$rank < ncol(design)
  n_pos_hc <- length(unique(d$donor_id[d$condition == "HC" & d$count > 0]))
  n_pos_ssc <- length(unique(d$donor_id[d$condition == "SSc" & d$count > 0]))
  out <- list(
    n_donors = length(unique(d$donor_id)),
    n_SSc = length(unique(d$donor_id[d$condition == "SSc"])),
    n_HC = length(unique(d$donor_id[d$condition == "HC"])),
    n_cohorts = nlevels(d$cohort),
    positive_donors_SSc = n_pos_ssc,
    positive_donors_HC = n_pos_hc,
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
  if (out$n_SSc < 2 || out$n_HC < 2) {
    out$model_message <- "fewer than two donors in one condition"
    return(out)
  }
  if (n_pos_hc < 2 || n_pos_ssc < 2) {
    out$model_message <- "not estimable: fewer than two positive donors in one condition"
    return(out)
  }
  if (rank_deficient) {
    out$model_message <- "rank-deficient condition-plus-cohort design"
    return(out)
  }
  fit <- tryCatch(
    glmmTMB(
      cbind(count, total - count) ~ condition + cohort,
      family = betabinomial(link = "logit"), data = d
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    out$model_message <- conditionMessage(fit)
    return(out)
  }
  coef_tab <- summary(fit)$coefficients$cond
  if (!"conditionSSc" %in% rownames(coef_tab)) {
    out$model_message <- "conditionSSc coefficient absent"
    return(out)
  }
  est <- coef_tab["conditionSSc", "Estimate"]
  se <- coef_tab["conditionSSc", "Std. Error"]
  p <- coef_tab["conditionSSc", "Pr(>|z|)"]
  pearson <- residuals(fit, type = "pearson")
  df_resid <- df.residual(fit)
  overdisp <- ifelse(is.finite(df_resid) && df_resid > 0, sum(pearson^2) / df_resid, NA)
  out$converged <- isTRUE(fit$sdr$pdHess) && fit$fit$convergence == 0
  out$pdHess <- isTRUE(fit$sdr$pdHess)
  out$dispersion <- unname(overdisp)
  out$estimate_log_or <- unname(est)
  out$se_log_or <- unname(se)
  out$OR <- unname(exp(est))
  out$CI_low <- unname(exp(est - 1.96 * se))
  out$CI_high <- unname(exp(est + 1.96 * se))
  out$p_value <- unname(p)
  out$model_message <- ifelse(out$converged, "ok", paste("convergence", fit$fit$convergence))
  out
}

fit_frame <- function(d, omitted_cohort = "none") {
  ans <- fit_one(d)
  data.frame(
    family = d$family[1], denominator = d$denominator[1], endpoint = d$endpoint[1],
    display = d$display[1], omitted_cohort = omitted_cohort,
    n_donors = ans$n_donors, n_SSc = ans$n_SSc, n_HC = ans$n_HC,
    n_cohorts = ans$n_cohorts,
    positive_donors_SSc = ans$positive_donors_SSc,
    positive_donors_HC = ans$positive_donors_HC,
    OR = ans$OR, CI_low = ans$CI_low, CI_high = ans$CI_high,
    p_value = ans$p_value, estimate_log_or = ans$estimate_log_or,
    se_log_or = ans$se_log_or, converged = ans$converged,
    pdHess = ans$pdHess, rank_deficient = ans$rank_deficient,
    dispersion_pearson = ans$dispersion, model_message = ans$model_message,
    stringsAsFactors = FALSE
  )
}

rows <- list()
keys <- unique(dat[, c("family", "denominator", "endpoint")])
for (i in seq_len(nrow(keys))) {
  key <- keys[i, ]
  d <- dat[
    dat$family == key$family & dat$denominator == key$denominator & dat$endpoint == key$endpoint,
  ]
  rows[[length(rows) + 1]] <- fit_frame(d)
}
res <- do.call(rbind, rows)
res$q_BH <- NA_real_
res$p_Holm <- NA_real_
families <- unique(res[, c("family", "denominator")])
for (i in seq_len(nrow(families))) {
  key <- families[i, ]
  ii <- which(
    res$family == key$family & res$denominator == key$denominator & is.finite(res$p_value)
  )
  if (length(ii)) {
    res$q_BH[ii] <- p.adjust(res$p_value[ii], method = "BH")
    res$p_Holm[ii] <- p.adjust(res$p_value[ii], method = "holm")
  }
}
write.csv(res, file.path(outdir, "Figure1c_glmmTMB_results.csv"), row.names = FALSE)

loco_rows <- list()
primary <- dat[dat$family == "fibro_primary" & dat$denominator == "all_skin_cells", ]
for (endpoint_name in unique(primary$endpoint)) {
  endpoint_data <- primary[primary$endpoint == endpoint_name, ]
  for (cohort_name in levels(droplevels(endpoint_data$cohort))) {
    d <- endpoint_data[endpoint_data$cohort != cohort_name, ]
    if (nrow(d)) {
      loco_rows[[length(loco_rows) + 1]] <- fit_frame(d, cohort_name)
    }
  }
}
if (length(loco_rows)) {
  write.csv(
    do.call(rbind, loco_rows), file.path(outdir, "Figure1c_leave_one_cohort_out.csv"),
    row.names = FALSE
  )
}

summary_keys <- c("family", "denominator", "endpoint", "display", "condition")
raw_summary <- aggregate(
  fraction ~ family + denominator + endpoint + display + condition,
  data = dat, FUN = median
)
names(raw_summary)[names(raw_summary) == "fraction"] <- "median"
for (stat_name in c("q1", "q3", "mean")) {
  stat_fun <- switch(
    stat_name,
    q1 = function(x) unname(quantile(x, 0.25)),
    q3 = function(x) unname(quantile(x, 0.75)),
    mean = mean
  )
  stat_table <- aggregate(
    fraction ~ family + denominator + endpoint + display + condition,
    data = dat, FUN = stat_fun
  )
  names(stat_table)[names(stat_table) == "fraction"] <- stat_name
  raw_summary <- merge(raw_summary, stat_table, by = summary_keys)
}
raw_n <- aggregate(
  donor_id ~ family + denominator + endpoint + display + condition,
  data = dat, FUN = function(x) length(unique(x))
)
names(raw_n)[names(raw_n) == "donor_id"] <- "n_donors"
raw_summary <- merge(
  raw_summary, raw_n,
  by = summary_keys
)
write.csv(raw_summary, file.path(outdir, "Figure1c_group_summary.csv"), row.names = FALSE)

cohort_summary <- aggregate(
  cbind(count, total) ~ family + denominator + endpoint + display + cohort + condition,
  data = dat, FUN = sum
)
cohort_summary$fraction <- cohort_summary$count / cohort_summary$total
write.csv(cohort_summary, file.path(outdir, "Figure1c_cohort_summary.csv"), row.names = FALSE)
