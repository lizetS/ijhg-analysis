# ===========================================================================
# 24_negative_binomial.R -- negative-binomial sensitivity for overdispersion.
#
# Reviewer 2, point 7. CARBayes has no negative-binomial likelihood, so this is
# a non-spatial negative-binomial GLM (MASS::glm.nb) with the same covariates
# and offset as the corresponding spatial model. It tests whether the fixed-
# effect estimates are sensitive to the Poisson variance assumption; it does
# not replace the spatial model.
#
# Output: output/tables/sensitivity_negative_binomial.csv
# ===========================================================================

RUN_SENSITIVITY <- TRUE
source("R/04_fit_models.R")

log_msg("--- 24_negative_binomial ---")

nb_one <- function(outcome, wave) {
  d <- DF %>% filter(covid_variant == wave) %>% arrange(Etiqueta)
  covs <- intersect(covariates_for(outcome, wave), names(d))
  f <- build_formula(covs, outcome)
  fit <- tryCatch(MASS::glm.nb(f, data = d),
                  error = function(e) {
                    log_msg("  glm.nb failed for ", outcome, "_", wave, ": ",
                            conditionMessage(e))
                    NULL
                  })
  if (is.null(fit)) return(NULL)
  b  <- stats::coef(fit)
  se <- sqrt(diag(stats::vcov(fit)))
  tibble(
    outcome = outcome, wave = WAVE_LABELS[[wave]],
    variable = names(b),
    rr = exp(b), lower = exp(b - 1.96 * se), upper = exp(b + 1.96 * se),
    theta = fit$theta, dispersion_glm = summary(fit)$dispersion,
    n = nrow(d)
  )
}

res <- bind_rows(lapply(c("incidence", "mortality"), function(oc) {
  waves <- if (oc == "incidence") WAVES else SEVERITY_WAVES
  bind_rows(lapply(waves, function(w) nb_one(oc, w)))
}))

write_csv_audit(res %>% filter(variable != "(Intercept)"),
                "sensitivity_negative_binomial.csv")
log_msg("24_negative_binomial complete: ", nrow(res), " rows")
