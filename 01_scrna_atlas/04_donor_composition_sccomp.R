suppressPackageStartupMessages({
  library(dplyr)
  library(sccomp)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: Rscript 04_donor_composition_sccomp.R counts.csv outdir")
}
counts_path <- args[[1]]
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# sccomp is the literature-anchored, sum-constrained beta-binomial sensitivity
# analysis (Mangiola et al., PNAS 2023). The frequentist glmmTMB model remains
# primary because it supplies cohort-adjusted odds ratios and Wald intervals.
dat <- read.csv(counts_path, stringsAsFactors = FALSE)
dat <- dat[dat$condition %in% c("HC", "SSc") & dat$total > 0, ]
dat$condition <- factor(dat$condition, levels = c("HC", "SSc"))
dat$cohort <- factor(dat$cohort)

fix_sccomp_cache <- function() {
  cache <- file.path(path.expand("~"), ".sccomp_models", as.character(packageVersion("sccomp")))
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  try(assignInNamespace("sccomp_stan_models_cache_dir", cache, ns = "sccomp"), silent = TRUE)
  cache
}

run_sccomp <- function(long, analysis_name) {
  long <- long %>%
    transmute(
      sample = donor_id,
      cell_group,
      count = as.integer(count),
      condition = droplevels(condition),
      cohort = droplevels(cohort)
    )
  if (any(long$count < 0)) stop(sprintf("negative counts in %s", analysis_name))
  totals <- aggregate(count ~ sample, data = long, FUN = sum)
  if (length(unique(totals$count)) < 2) {
    message(sprintf("%s: all donors have the same total", analysis_name))
  }
  draws_dir <- file.path(outdir, paste0("sccomp_draws_", analysis_name))
  fit <- sccomp_estimate(
    long,
    formula_composition = ~ condition + cohort,
    formula_variability = ~ condition,
    .sample = sample,
    .cell_group = cell_group,
    .abundance = count,
    cores = 4,
    inference_method = "pathfinder",
    percent_false_positive = 5,
    output_directory = draws_dir,
    mcmc_seed = 20260723,
    verbose = FALSE
  )
  tested <- sccomp_test(
    fit,
    percent_false_positive = 5,
    test_composition_above_logit_fold_change = 0
  )
  condition_rows <- tested %>%
    filter(!is.na(factor) & grepl("^condition", factor)) %>%
    mutate(analysis = analysis_name, .before = 1) %>%
    select(any_of(c(
      "analysis", "cell_group", "factor",
      "c_effect", "c_lower", "c_upper", "c_pH0", "c_FDR",
      "v_effect", "v_lower", "v_upper", "v_pH0", "v_FDR"
    )))
  if (!nrow(condition_rows)) {
    stop(sprintf("condition coefficient absent from sccomp result for %s", analysis_name))
  }
  condition_rows
}

fix_sccomp_cache()

# Major cell types are mutually exclusive and sum to every retained skin cell.
major <- dat %>%
  filter(family == "major_celltype", denominator == "all_skin_cells") %>%
  transmute(donor_id, cohort, condition, cell_group = endpoint, count)
major_result <- run_sccomp(major, "major_celltypes")

# The seven fibroblast states are made sum-constrained against all other skin
# cells. The total Fibroblast endpoint is omitted because it overlaps its states.
primary <- dat %>%
  filter(family == "fibro_primary", denominator == "all_skin_cells")
subtypes <- c(
  "Myofibroblast", "SFRP4_proFib", "SFRP2_DPP4", "Adipogenic",
  "FMO1_LSP1", "LGR5_Gur", "Inflammatory"
)
state_rows <- primary %>%
  filter(endpoint %in% subtypes) %>%
  transmute(donor_id, cohort, condition, cell_group = endpoint, count, total)
other_rows <- state_rows %>%
  group_by(donor_id, cohort, condition) %>%
  summarise(
    count = first(total) - sum(count),
    .groups = "drop"
  ) %>%
  mutate(cell_group = "Other_skin_cells")
state_long <- bind_rows(
  state_rows %>% select(donor_id, cohort, condition, cell_group, count),
  other_rows %>% select(donor_id, cohort, condition, cell_group, count)
)
state_result <- run_sccomp(state_long, "fibroblast_state_contributions")

result <- bind_rows(major_result, state_result)
write.csv(result, file.path(outdir, "Figure1c_sccomp_sensitivity.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(outdir, "Figure1c_sccomp_sessionInfo.txt"))
message(sprintf("sccomp sensitivity complete: %s", outdir))
