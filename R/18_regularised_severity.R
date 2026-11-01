# ===========================================================================
# 18_regularised_severity.R -- sensitivity analysis for the severity models.
#
# Why: severity models carry 18 covariates (12 structural + 6 comorbidities).
# Early (48 deaths) and Wuhan (20 deaths) have too few events for that. No
# covariate is dropped (a priori or after looking at results). Instead:
#   A. collinearity is REPORTED (VIF of all 18 covariates; no variable removed)
#   B. every severity model is refitted on the z-score scale with a shrinkage
#      prior on the regression coefficients, N(0, 1) on the log-RR per SD
#      (95% prior range of the RR per SD: 0.14 to 7.1), and compared with the
#      same model with a vague prior N(0, 1e5).
# The primary models (min-max scale, CARBayes default priors) are unchanged.
#
# Outputs: output/tables/collinearity_vif_severity.csv
#          output/tables/sensitivity_regularised_severity.csv
#          output/tables/sensitivity_regularised_summary.csv
#          output/tables/sensitivity_regularised_diagnostics.csv
# ===========================================================================

source("config.R")
source(file.path(ROOT, "R", "00_setup.R"))

RUN_SENSITIVITY <- TRUE
source(file.path(ROOT, "R", "04_fit_models.R"))   # fit_model(), DF, W_QUEEN

log_msg("--- 18_regularised_severity ---")
label_of <- function(v) if (v %in% names(LABELS)) unname(LABELS[v]) else v
CFG$thin_severity <- 5L

# ---- A. VIF, reported without removing anything -----------------------------
vif_tab <- function(d, vars) {
  bind_rows(lapply(vars, function(v) {
    f <- stats::as.formula(paste(v, "~", paste(setdiff(vars, v), collapse = "+")))
    r2 <- summary(stats::lm(f, data = d))$r.squared
    tibble(variable = v, vif = ifelse(r2 >= 1, Inf, 1 / (1 - r2)))
  }))
}
d1 <- DF %>% filter(covid_variant == WAVES[1])
vif_sets <- list(
  structural_12 = intersect(COV_INCIDENCE, names(d1)),
  severity_18   = intersect(COV_SEVERITY,  names(d1))
)
vif_all <- bind_rows(lapply(names(vif_sets), function(s)
  vif_tab(d1, vif_sets[[s]]) %>% mutate(set = s, .before = 1))) %>%
  mutate(label = vapply(variable, label_of, character(1)),
         above_5 = vif > 5, above_10 = vif > 10)
# mortality models also carry the wave-specific incidence covariate
vif_mort <- bind_rows(lapply(SEVERITY_WAVES, function(w) {
  dw <- DF %>% filter(covid_variant == w)
  vif_tab(dw, intersect(c(COV_SEVERITY, COV_WAVE_VARYING), names(dw))) %>%
    mutate(set = paste0("mortality_with_incidence_", w), .before = 1,
           label = vapply(variable, label_of, character(1)),
           above_5 = vif > 5, above_10 = vif > 10)
}))
vif_all <- bind_rows(vif_all, vif_mort)
write_csv_audit(vif_all, "collinearity_vif_severity.csv")
log_msg("max VIF structural_12 = ", round(max(vif_all$vif[vif_all$set == "structural_12"]), 2),
        "; severity_18 = ", round(max(vif_all$vif[vif_all$set == "severity_18"]), 2))

# ---- B. z-scores + shrinkage prior ------------------------------------------
raw <- read_rds("analysis_dataset.rds") %>%
  st_drop_geometry() %>% arrange(Etiqueta, covid_variant)
covs_z <- setdiff(intersect(COV_SEVERITY, names(raw)), COV_WAVE_VARYING)
for (v in covs_z) {
  x <- raw[[v]][raw$covid_variant == WAVES[1]]          # 82 area values
  raw[[v]] <- (raw[[v]] - mean(x)) / sd(x)
}
if ("incidence_cov" %in% names(raw)) {      # wave-specific covariate: z-score within each wave
  raw <- raw %>% group_by(covid_variant) %>%
    mutate(incidence_cov = (incidence_cov - mean(incidence_cov)) / sd(incidence_cov)) %>%
    ungroup()
}
DF_primary <- DF
DF <- raw                                               # fit_model() reads DF

prior_sets <- list(
  vague_z  = function(p) list(prior.var.beta = c(1e5, rep(1e5, p - 1))),
  shrink_z = function(p) list(prior.var.beta = c(1e5, rep(1, p - 1)))   # N(0,1) on slopes
)

res <- list(); dg <- list()
for (outcome in c("mortality", "cfr")) {
  for (w in SEVERITY_WAVES) {
    covs <- covariates_for(outcome, w)
    nd <- sum(DF$deaths_per_variant[DF$covid_variant == w])
    for (pn in names(prior_sets)) {
      tag <- paste(outcome, w, pn, sep = "_")
      log_msg("fitting ", tag, " (", nd, " deaths, ", length(covs), " covariates)")
      r <- tryCatch(fit_model(outcome, w, covs, priors = prior_sets[[pn]](length(covs) + 1)),
                    error = function(e) { log_msg("  FAILED: ", conditionMessage(e)); NULL })
      if (is.null(r)) next
      res[[tag]] <- r$coef %>% mutate(prior = pn, n_deaths = nd,
                                      events_per_variable = nd / length(covs))
      dg[[tag]]  <- r$diag %>% mutate(prior = pn, n_deaths = nd)
      rm(r); invisible(gc())
    }
  }
}
DF <- DF_primary

coef_z <- bind_rows(res) %>%
  filter(variable != "(Intercept)") %>%
  mutate(excludes_1 = lower > 1 | upper < 1,
         ci_width_log = log(upper) - log(lower))

# significance in the primary (min-max, default prior) models, for comparison
prim <- utils::read.csv(file.path(DIR$tables, "coefficients.csv"), stringsAsFactors = FALSE) %>%
  filter(outcome %in% c("mortality", "cfr")) %>%
  transmute(outcome, wave, variable, rr_primary = rr,
            excludes_1_primary = lower > 1 | upper < 1)
coef_z <- coef_z %>% left_join(prim, by = c("outcome", "wave", "variable"))
write_csv_audit(coef_z, "sensitivity_regularised_severity.csv")
write_csv_audit(bind_rows(dg), "sensitivity_regularised_diagnostics.csv")

summ <- coef_z %>%
  group_by(outcome, wave, n_deaths, events_per_variable, prior) %>%
  summarise(n_coef = n(), n_excl_1 = sum(excludes_1),
            median_ci_width_log = median(ci_width_log), .groups = "drop") %>%
  left_join(coef_z %>% filter(prior == "vague_z") %>%
              group_by(outcome, wave) %>%
              summarise(n_excl_1_primary = sum(excludes_1_primary, na.rm = TRUE), .groups = "drop"),
            by = c("outcome", "wave"))
write_csv_audit(summ, "sensitivity_regularised_summary.csv")
print(summ)
log_msg("18_regularised_severity complete")
