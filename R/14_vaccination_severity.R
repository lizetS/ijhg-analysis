# ===========================================================================
# 14_vaccination_severity.R -- sensitivity analysis: add area-level vaccination
# coverage to the Delta mortality and case-fatality models.
#
# Reviewer 2 (point 14) asks why vaccination was modelled for incidence but not
# for severity. Coverage exists only in the Delta wave (`vacunados` is empty in
# every other wave), and Delta has 524 deaths, so the model is estimable. This
# answers the question with data rather than with a rationale.
#
# Output: output/tables/sensitivity_vaccination_severity.csv
# ===========================================================================

source("config.R")
source(file.path(ROOT, "R", "00_setup.R"))

RUN_SENSITIVITY <- TRUE
source(file.path(ROOT, "R", "04_fit_models.R"))

log_msg("--- 14_vaccination_severity ---")
CFG$thin_severity <- 5L

out <- list()
for (outcome in c("mortality", "cfr")) {
  covs <- c(covariates_for(outcome, "Delta"), COV_VACCINE)
  log_msg("fitting ", outcome, "_Delta with vaccination coverage (",
          length(covs), " covariates)")
  r <- tryCatch(fit_model(outcome, "Delta", covs),
                error = function(e) { log_msg("  FAILED: ", conditionMessage(e)); NULL })
  if (is.null(r)) next
  log_msg("  rhat ", round(r$diag$rhat_max, 4), "; ess ", round(r$diag$ess_min))
  out[[outcome]] <- r$coef %>% filter(variable %in% c("cobertura_vac", "(Intercept)"))
  save_rds(r$obj, paste0(outcome, "_Delta_vaccination.rds"), DIR$models)
}

res <- bind_rows(out)
write_csv_audit(res, "sensitivity_vaccination_severity.csv")
cat("\nCoeficiente de cobertura vacunal en los modelos de gravedad de Delta:\n")
print(as.data.frame(res %>% filter(variable == "cobertura_vac")),
      row.names = FALSE)
