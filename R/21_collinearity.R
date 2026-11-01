# ===========================================================================
# 21_collinearity.R -- collinearity diagnostics for every model (Reviewer 2,
# point 2.5).
#
# For each outcome and wave, on the covariate set that model actually uses:
#   - variance inflation factor of each covariate (1 / (1 - R2) of the
#     auxiliary regression on the remaining covariates), and
#   - its strongest absolute pairwise correlation, and with which covariate.
#
# No model has to be fitted: these are linear regressions on the covariates.
#
# Outputs: output/tables/collinearity.csv          one row per covariate
#          output/tables/collinearity_summary.csv one row per model
# ===========================================================================

source("config.R")
source(file.path(ROOT, "R", "00_setup.R"))

log_msg("--- 21_collinearity ---")
dat <- read_rds("analysis_standardised.rds") %>% st_drop_geometry()

vif_one <- function(d, vars) {
  vapply(vars, function(v) {
    f <- stats::as.formula(paste(v, "~", paste(setdiff(vars, v), collapse = " + ")))
    r2 <- summary(stats::lm(f, data = d))$r.squared
    if (r2 >= 1) Inf else 1 / (1 - r2)
  }, numeric(1))
}

rows <- list(); summary_rows <- list()

for (oc in OUTCOMES) {
  for (w in WAVES) {
    d <- dat %>% filter(covid_variant == w)
    vars <- intersect(covariates_for(oc, w), names(d))
    if (length(vars) < 2) next

    v <- vif_one(d, vars)
    M <- stats::cor(d[, vars])
    diag(M) <- NA

    A <- abs(M)
    best_r <- numeric(length(vars))
    best_with <- character(length(vars))
    for (k in seq_along(vars)) {
      row <- A[vars[k], ]
      row <- row[!is.na(row)]
      best_r[k] <- if (length(row)) max(row) else NA_real_
      best_with[k] <- if (length(row)) names(row)[which.max(row)] else NA_character_
    }

    rows[[length(rows) + 1]] <- tibble(
      outcome = oc, wave = WAVE_LABELS[[w]], n_covariates = length(vars),
      variable = vars,
      vif = round(unname(v[vars]), 3),
      max_abs_r = round(best_r, 3),
      max_abs_r_with = best_with
    )

    summary_rows[[length(summary_rows) + 1]] <- tibble(
      outcome = oc, wave = WAVE_LABELS[[w]], n_covariates = length(vars),
      vif_max = round(max(v), 3),
      vif_max_variable = vars[which.max(v)],
      n_vif_gt_5 = sum(v > 5), n_vif_gt_10 = sum(v > 10),
      max_abs_r = round(max(abs(M), na.rm = TRUE), 3)
    )
  }
}

written <- bind_rows(rows)
summary_tab <- bind_rows(summary_rows)

write_csv_audit(written, "collinearity.csv")
write_csv_audit(summary_tab, "collinearity_summary.csv")

log_msg("21_collinearity complete: VIF maximo ",
        round(max(summary_tab$vif_max), 2), " (",
        summary_tab$vif_max_variable[which.max(summary_tab$vif_max)], "); ",
        sum(summary_tab$n_vif_gt_5 > 0), " modelos con algun VIF > 5")
print(summary_tab %>% select(outcome, wave, n_covariates, vif_max, vif_max_variable, n_vif_gt_5))
