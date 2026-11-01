# ===========================================================================
# 25_extra_outputs.R -- items requested in review that are derived from the
# stored models and the analysis dataset (no refitting).
#
#   coefficients_99ci.csv            99% credible intervals (R1 1.23)
#   morans_I_alternative_weights.csv Moran's I under four neighbourhoods (R2 2.9)
#   vaccination_coverage_summary.csv coverage definition and spatial variation (R2 2.14)
#   Table3_waves.csv / .xlsx         wave periods with dates (R1 1.5)
#   TableS5_diagnostics.xlsx         convergence (R1 1.14, R2 2.6)
#   TableS6_vaccination_severity.xlsx
#   TableS7_model_validation.xlsx    dispersion, PPC, residual Moran (R2 2.6, 2.7)
# ===========================================================================

source("config.R")
source(file.path(ROOT, "R", "00_setup.R"))

log_msg("--- 25_extra_outputs ---")
dat <- read_rds("analysis_dataset.rds")

# ---- 1. 99% credible intervals from the stored traces ----------------------
files <- list.files(DIR$models, "\\.rds$", full.names = TRUE)
ci99 <- bind_rows(lapply(files, function(f) {
  o <- readRDS(f)
  if (is.null(o$trace)) return(NULL)
  tr <- as.matrix(o$trace[, setdiff(colnames(o$trace), "chain"), drop = FALSE])
  tibble(outcome = o$outcome, wave = WAVE_LABELS[[o$wave]],
         variable = colnames(tr),
         lower99 = exp(apply(tr, 2, stats::quantile, 0.005)),
         upper99 = exp(apply(tr, 2, stats::quantile, 0.995)),
         n_draws = nrow(tr))
}))
write_csv_audit(ci99, "coefficients_99ci.csv")

# ---- 2. Moran's I under alternative neighbourhood definitions -------------
SP <- dat %>% filter(covid_variant == WAVES[1]) %>% arrange(Etiqueta)
coords <- sf::st_coordinates(sf::st_centroid(SP))
NB <- spdep::poly2nb(SP, queen = TRUE)
W <- list(
  queen  = spdep::nb2mat(NB, style = "B"),
  rook   = spdep::nb2mat(spdep::poly2nb(SP, queen = FALSE), style = "B"),
  queen2 = spdep::nb2mat(spdep::nblag_cumul(spdep::nblag(NB, 2L)), style = "B",
                         zero.policy = TRUE),
  knn4   = spdep::nb2mat(spdep::knn2nb(spdep::knearneigh(coords, k = 4)), style = "B")
)
W$knn4 <- 1L * ((W$knn4 + t(W$knn4)) > 0)

moran_tab <- bind_rows(lapply(names(W), function(nm) {
  lw <- spdep::mat2listw(W[[nm]], style = "W")
  bind_rows(lapply(OUTCOMES, function(oc) {
    var <- switch(oc, incidence = "incidence_rate", mortality = "mortality_rate",
                  cfr = "lethality")
    bind_rows(lapply(WAVES, function(w) {
      d <- dat %>% filter(covid_variant == w) %>% arrange(Etiqueta)
      mc <- spdep::moran.mc(d[[var]], lw, nsim = 1001, zero.policy = TRUE)
      tibble(weights = nm, outcome = oc, wave = WAVE_LABELS[[w]],
             moran_I = unname(mc$statistic), p_value = mc$p.value)
    }))
  }))
}))
write_csv_audit(moran_tab, "morans_I_alternative_weights.csv")

# ---- 3. vaccination coverage: definition and spatial variation ------------
vacc <- dat %>% st_drop_geometry() %>%
  filter(!is.na(cobertura_vac)) %>%
  group_by(covid_variant) %>%
  summarise(n_areas = n(), mean_pct = mean(cobertura_vac), sd_pct = sd(cobertura_vac),
            min_pct = min(cobertura_vac), q25 = quantile(cobertura_vac, 0.25),
            median_pct = median(cobertura_vac), q75 = quantile(cobertura_vac, 0.75),
            max_pct = max(cobertura_vac), .groups = "drop") %>%
  mutate(wave = WAVE_LABELS[covid_variant], .after = covid_variant) %>%
  select(-covid_variant)
write_csv_audit(vacc, "vaccination_coverage_summary.csv")

# ---- 4. Table 3 with wave dates -------------------------------------------
wave_tab <- dat %>% st_drop_geometry() %>%
  group_by(wave = covid_variant) %>%
  summarise(start = min(as.Date(start_date), na.rm = TRUE), n = n(),
            pop = sum(poblacion), cases = sum(cases_per_variant),
            deaths = sum(deaths_per_variant), .groups = "drop") %>%
  mutate(wave = factor(wave, levels = WAVES)) %>% arrange(wave) %>%
  mutate(end = as.Date(dplyr::lead(start) - 1),
         end = as.Date(ifelse(is.na(end), as.Date("2022-05-31"), end),
                       origin = "1970-01-01"),
         label = WAVE_LABELS[as.character(wave)],
         lineage = c("Before variant characterisation", "Ancestral (Wuhan) lineage",
                     "D614G", "B.1.351", "B.1.617.2", "BA.1", "BA.2"),
         cum_incidence_per_100k = round(cases / pop * 1e5, 1),
         cfr_pct = round(deaths / cases * 100, 2)) %>%
  select(label, lineage, start, end, n, cases, deaths,
         cum_incidence_per_100k, cfr_pct)
write_csv_audit(wave_tab, "Table3_waves.csv")
writexl::write_xlsx(wave_tab, file.path(DIR$tables, "Table3_waves.xlsx"))

# ---- 5-7. formatted supplementary tables ----------------------------------
diag_tab <- utils::read.csv(file.path(DIR$tables, "diagnostics.csv"),
                            stringsAsFactors = FALSE)
s5 <- diag_tab[, intersect(c("outcome", "wave", "n_chains", "n_iter", "burnin", "thin",
                             "sampler", "rhat_max", "ess_min", "rho_mean", "tau2_mean",
                             "status"), names(diag_tab))]
writexl::write_xlsx(s5, file.path(DIR$tables, "TableS5_diagnostics.xlsx"))

f_s6 <- file.path(DIR$tables, "sensitivity_vaccination_severity.csv")
if (file.exists(f_s6))
  writexl::write_xlsx(utils::read.csv(f_s6, stringsAsFactors = FALSE),
                      file.path(DIR$tables, "TableS6_vaccination_severity.xlsx"))

s7 <- diag_tab[, intersect(c("outcome", "wave", "pearson_dispersion", "pearson_ppp",
                             "moran_residual", "moran_p", "sparse_outcome", "status"),
                           names(diag_tab))]
writexl::write_xlsx(s7, file.path(DIR$tables, "TableS7_model_validation.xlsx"))

log_msg("25_extra_outputs complete")
