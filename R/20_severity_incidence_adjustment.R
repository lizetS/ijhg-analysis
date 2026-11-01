# ===========================================================================
# 20_severity_incidence_adjustment.R -- sensitivity: adjustment for incidence
# in the severity models.
#
# Primary models (config.R): mortality adjusted for the incidence of the same
# wave, case fatality not adjusted. This script fits all four combinations for
# the five severity waves and compares the coefficients:
#   mortality without incidence  (shows what the adjustment changes)
#   mortality with incidence     (primary)
#   case fatality without incidence (primary)
#   case fatality with incidence (sensitivity; note mathematical coupling,
#                                 cases are both the offset and the base of
#                                 the incidence covariate)
# Output: output/tables/sensitivity_incidence_adjustment.csv
# ===========================================================================

source("config.R")
source(file.path(ROOT, "R", "00_setup.R"))
RUN_SENSITIVITY <- TRUE
source(file.path(ROOT, "R", "04_fit_models.R"))   # fit_model(), DF
log_msg("--- 20_severity_incidence_adjustment ---")
CFG$thin_severity <- 5L

specs <- list(
  list(outcome = "mortality", model = "without_incidence", covs = COV_SEVERITY),
  list(outcome = "mortality", model = "with_incidence",    covs = c(COV_SEVERITY, COV_WAVE_VARYING)),
  list(outcome = "cfr",       model = "without_incidence", covs = COV_SEVERITY),
  list(outcome = "cfr",       model = "with_incidence",    covs = c(COV_SEVERITY, COV_WAVE_VARYING))
)

res <- list()
for (sp in specs) {
  for (w in SEVERITY_WAVES) {
    tag <- paste(sp$outcome, w, sp$model, sep = "_")
    log_msg("fitting ", tag)
    r <- tryCatch(fit_model(sp$outcome, w, sp$covs),
                  error = function(e) { log_msg("  FAILED: ", conditionMessage(e)); NULL })
    if (is.null(r)) next
    res[[tag]] <- r$coef %>% mutate(model = sp$model, .before = 1)
    rm(r); invisible(gc())
  }
}

long <- bind_rows(res) %>% filter(variable != "(Intercept)") %>%
  mutate(excludes_1 = lower > 1 | upper < 1)
wide <- long %>%
  select(outcome, wave, variable, model, rr, lower, upper, excludes_1) %>%
  tidyr::pivot_wider(names_from = model, values_from = c(rr, lower, upper, excludes_1)) %>%
  mutate(status_changes = excludes_1_without_incidence != excludes_1_with_incidence)
write_csv_audit(wide, "sensitivity_incidence_adjustment.csv")

print(wide %>% filter(variable == "incidence_cov") %>%
        select(outcome, wave, rr_with_incidence, lower_with_incidence, upper_with_incidence))
log_msg(sum(wide$status_changes, na.rm = TRUE), " coefficients change null-interval status with the adjustment")
log_msg("20_severity_incidence_adjustment complete")
