# ===========================================================================
# 08_sensitivity.R -- the analyses requested by the reviewers.
#
#   A. Incidence models with chronic-disease prevalence   (R1 1.11, R2 5) -> Table S4
#   B. Case-fatality models with a binomial likelihood    (R2 3, 7)
#   C. Alternative spatial weights: rook, 2nd-order queen, k-NN (R2 9)
#   D. Alternative priors for the CAR variance             (R2 8)
#   E. Collinearity: correlations and variance inflation   (R2 5)
#   F. Wave-definition robustness                          (R2 15)
#
# This script is expensive: A and B refit 7 models each, C and D refit the
# incidence models once per alternative. Run it deliberately, not as part of
# every pipeline run.
#
# Outputs (in output/tables/): TableS4_incidence_comorbid.csv,
#   sensitivity_binomial_cfr.csv, sensitivity_weights.csv,
#   sensitivity_priors.csv, collinearity.csv, wave_robustness.csv
# ===========================================================================


source("config.R")
source(file.path(ROOT, "R", "00_setup.R"))

RUN_SENSITIVITY <- TRUE
source(file.path(ROOT, "R", "04_fit_models.R"))   # defines fit_model only

log_msg("--- 08_sensitivity ---")
prim <- utils::read.csv(file.path(DIR$tables, "coefficients.csv"),
                        stringsAsFactors = FALSE) %>%
  filter(variable != "(Intercept)")

# ============================ A. comorbidities =============================
log_msg("A. incidence models with chronic-disease prevalence")
comorb <- bind_rows(lapply(WAVES, function(w) {
  covs <- c(covariates_for("incidence", w), COV_COMORBID)
  covs <- intersect(covs, names(DF))
  r <- tryCatch(fit_model("incidence", w, covs)$coef,
                error = function(e) NULL)
  if (is.null(r)) { log_msg("  failed: ", w); return(NULL) }
  r
}))
write_csv_audit(comorb, "TableS4_incidence_comorbid.csv")

if (nrow(comorb)) {
  cmp <- comorb %>%
    filter(variable %in% COV_INCIDENCE) %>%
    select(outcome, wave, variable, rr_adj = rr, lo_adj = lower, hi_adj = upper) %>%
    inner_join(prim %>% filter(outcome == "incidence") %>%
                 select(wave, variable, rr_prim = rr, lo_prim = lower, hi_prim = upper),
               by = c("wave", "variable")) %>%
    mutate(
      direction_unchanged = sign(log(rr_adj)) == sign(log(rr_prim)),
      still_excludes_1    = lo_adj > 1 | hi_adj < 1,
      log_rr_change       = round(log(rr_adj) - log(rr_prim), 3)
    ) %>%
    mutate(variable = vapply(variable, label_of, character(1)))
  write_csv_audit(cmp, "TableS4_comparison.csv")
  log_msg("  comparison written; associations losing exclusion of 1: ",
          sum(!cmp$still_excludes_1 & (cmp$lo_prim > 1 | cmp$hi_prim < 1)))
}

# ====================== A2. reduced covariate set ==========================
# The 9-covariate structural set printed in Table 1 of the submitted
# manuscript, refitted as a sensitivity analysis (Reviewer 2, point 5).
log_msg("A2. incidence models with the reduced covariate set")
reduced <- bind_rows(lapply(WAVES, function(w) {
  covs <- COV_INCIDENCE_REDUCED
  if (w == "Delta") covs <- c(covs, COV_VACCINE)
  r <- tryCatch(fit_model("incidence", w, intersect(covs, names(DF)))$coef,
                error = function(e) NULL)
  if (is.null(r)) return(NULL)
  r %>% mutate(model = "reduced")
}))
if (nrow(reduced)) {
  full <- prim %>% filter(outcome == "incidence") %>% mutate(model = "full")
  cmp <- bind_rows(full, reduced) %>%
    select(model, wave, variable, rr, lower, upper) %>%
    tidyr::pivot_wider(names_from = model, values_from = c(rr, lower, upper))
  write_csv_audit(cmp, "sensitivity_reduced_covariates.csv")
}

# ============================ B. binomial CFR ==============================
log_msg("B. case-fatality models with a binomial likelihood")
binom <- bind_rows(lapply(SEVERITY_WAVES, function(w) {
  r <- tryCatch(fit_model("cfr", w, covariates_for("cfr", w), family = "binomial")$coef,
                error = function(e) { log_msg("  failed: ", w, " ", conditionMessage(e)); NULL })
  if (is.null(r)) return(NULL)
  r %>% mutate(model = "binomial")
}))
if (nrow(binom)) {
  pois <- prim %>% filter(outcome == "cfr") %>% mutate(model = "poisson")
  cmp <- bind_rows(pois, binom) %>%
    select(model, wave, variable, rr, lower, upper) %>%
    tidyr::pivot_wider(names_from = model,
                       values_from = c(rr, lower, upper)) %>%
    # the binomial model has a logit link: its coefficients are ODDS ratios, the
    # primary Poisson model (offset log cases) gives case-fatality RATIOS. Compare
    # direction and null-interval status; magnitudes differ by construction.
    rename(or_binomial = rr_binomial, lower_or_binomial = lower_binomial,
           upper_or_binomial = upper_binomial) %>%
    mutate(same_direction = (rr_poisson >= 1) == (or_binomial >= 1),
           excludes_1_poisson = lower_poisson > 1 | upper_poisson < 1,
           excludes_1_binomial = lower_or_binomial > 1 | upper_or_binomial < 1,
           same_status = excludes_1_poisson == excludes_1_binomial)
  write_csv_audit(cmp, "sensitivity_binomial_cfr.csv")
}

# ============================ C. spatial weights ===========================
log_msg("C. alternative spatial weights")
coords <- st_coordinates(st_centroid(SP))
W_alt <- list(
  queen     = W_QUEEN,
  rook      = nb2mat(poly2nb(SP, queen = FALSE), style = "B"),
  queen2    = nb2mat(nblag_cumul(nblag(NB, 2L)), style = "B", zero.policy = TRUE),
  knn4      = nb2mat(knn2nb(knearneigh(coords, k = 4)), style = "B")
)
# CARBayes requires a symmetric W; the k-NN graph is symmetrised.
W_alt <- lapply(W_alt, function(M) { M <- 1L * ((M + t(M)) > 0); M })  # symmetrise k-NN

weights_res <- bind_rows(lapply(names(W_alt), function(nm) {
  Wm <- W_alt[[nm]]
  if (any(rowSums(Wm) == 0)) { log_msg("  ", nm, ": isolated areas, skipped"); return(NULL) }
  bind_rows(lapply(WAVES, function(w) {
    r <- tryCatch(fit_model("incidence", w, covariates_for("incidence", w), W = Wm)$coef,
                  error = function(e) NULL)
    if (is.null(r)) return(NULL)
    r %>% mutate(weights = nm)
  }))
}))
if (nrow(weights_res))
  write_csv_audit(weights_res %>% filter(variable != "(Intercept)"),
                  "sensitivity_weights.csv")

# ============================ D. priors ====================================
log_msg("D. alternative priors for the CAR variance")
prior_sets <- list(
  # CARBayes defaults = what the primary models use: beta ~ N(0, 1e5) (variance),
  # tau2 ~ Inverse-Gamma(1, 0.01), rho ~ Uniform(0, 1)
  default  = list(),
  ig_1_001 = list(prior.tau2 = c(1, 0.01)),
  ig_half  = list(prior.tau2 = c(0.5, 0.0005)),
  # priors stated in the earlier README / draft Methods: N(0, 100) on beta
  # (variance 100) and Gamma(0.5, 0.0005) on the precision of tau2
  declared_N100_G05 = function(p) list(prior.var.beta = rep(100, p),
                                       prior.tau2 = c(0.5, 0.0005)),
  beta_N100 = function(p) list(prior.var.beta = rep(100, p))
)
prior_res <- bind_rows(lapply(names(prior_sets), function(nm) {
  bind_rows(lapply(WAVES, function(w) {
    covs <- covariates_for("incidence", w)
    pr <- prior_sets[[nm]]
    if (is.function(pr)) pr <- pr(length(covs) + 1)
    r <- tryCatch(fit_model("incidence", w, covs, priors = pr)$coef,
                  error = function(e) NULL)
    if (is.null(r)) return(NULL)
    r %>% mutate(prior = nm)
  }))
}))
if (nrow(prior_res))
  write_csv_audit(prior_res %>% filter(variable != "(Intercept)"),
                  "sensitivity_priors.csv")

# ============================ E. collinearity ==============================
log_msg("E. collinearity diagnostics")
vif_tab <- function(d, vars) {
  out <- lapply(vars, function(v) {
    f <- stats::as.formula(paste(v, "~", paste(setdiff(vars, v), collapse = "+")))
    fit <- stats::lm(f, data = d)
    r2 <- summary(fit)$r.squared
    tibble(variable = v, r2 = r2, vif = ifelse(r2 >= 1, Inf, 1 / (1 - r2)))
  })
  bind_rows(out)
}
col_tab <- bind_rows(lapply(WAVES, function(w) {
  d <- DF %>% filter(covid_variant == w)
  covs <- intersect(covariates_for("incidence", w), names(d))
  vif_tab(d, covs) %>% mutate(wave = WAVE_LABELS[[w]], .before = 1)
})) %>% mutate(variable = vapply(variable, label_of, character(1)))
write_csv_audit(col_tab, "collinearity.csv")
print(col_tab %>% group_by(wave) %>% summarise(max_vif = max(vif, na.rm = TRUE)))

# ============================ F. wave robustness ===========================
log_msg("F. wave-definition robustness")
linelist <- suppressMessages(readxl::read_excel(FILES$daily, .name_repair = "minimal"))
daily <- linelist %>%
  filter(!is.na(date_conf)) %>%
  mutate(date = as.Date(date_conf)) %>%
  count(date, name = "cases")

wave_algo <- function(daily, window = 15, interval = 14, r2_min = 0.75) {
  d <- daily %>% arrange(date) %>%
    tidyr::complete(date = seq(min(date), max(date), by = "day"),
                    fill = list(cases = 0))
  d$ma <- as.numeric(stats::filter(d$cases, rep(1 / window, window), sides = 1))
  d <- d %>% filter(!is.na(ma))
  d$seg <- ceiling(seq_len(nrow(d)) / interval)
  sym <- vapply(split(seq_len(nrow(d)), d$seg), function(i) {
    if (length(i) < 5 || all(is.na(d$ma[i]))) return("S")
    fit <- stats::lm(d$ma[i] ~ seq_along(i))
    sl <- coef(fit)[2]; r2 <- summary(fit)$r.squared
    if (sl > 0 && r2 >= r2_min) "R" else if (sl < 0 && r2 >= r2_min) "F" else "S"
  }, character(1))
  s <- paste(sym, collapse = "")
  m <- gregexpr("R+S*F+", s, perl = TRUE)[[1]]
  n <- if (m[1] == -1) 0 else length(m)
  tibble(window = window, interval = interval, r2_min = r2_min, n_waves = n)
}

grid <- expand.grid(window = c(7, 15, 21), interval = c(7, 14, 21),
                    r2_min = c(0.60, 0.75, 0.90))
wave_rob <- bind_rows(lapply(seq_len(nrow(grid)), function(i)
  wave_algo(daily, grid$window[i], grid$interval[i], grid$r2_min[i])))
write_csv_audit(wave_rob, "wave_robustness.csv")
print(wave_rob %>% count(n_waves))

log_msg("08_sensitivity complete")
