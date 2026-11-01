# ===========================================================================
# 02_waves_descriptives.R -- wave table (manuscript Table 3) and Moran's I.
# Outputs: output/tables/table3_waves.csv
#          output/tables/morans_I.csv
# ===========================================================================


source("config.R")
source(file.path(ROOT, "R", "00_setup.R"))

log_msg("--- 02_waves_descriptives ---")
dat <- read_rds("analysis_dataset.rds")

# ---- wave periods ----------------------------------------------------------
if ("start_date" %in% names(dat) && any(!is.na(dat$start_date))) {
  wave_tab <- dat %>%
    st_drop_geometry() %>%
    group_by(wave = covid_variant) %>%
    summarise(
      start = min(as.Date(start_date), na.rm = TRUE),
      n     = n(),
      pop   = sum(poblacion),
      cases = sum(cases_per_variant),
      deaths = sum(deaths_per_variant),
      .groups = "drop"
    )
} else {
  wave_tab <- dat %>%
    st_drop_geometry() %>%
    group_by(wave = covid_variant) %>%
    summarise(n = n(), pop = sum(poblacion), cases = sum(cases_per_variant),
              deaths = sum(deaths_per_variant), .groups = "drop")
}

lineage <- c(
  Wuhan1 = "Before variant characterisation", Wuhan2 = "Ancestral (Wuhan) lineage",
  D614G = "D614G", Beta = "B.1.351", Delta = "B.1.617.2",
  Omicron1 = "BA.1", Omicron2 = "BA.2"
)

wave_tab <- wave_tab %>%
  mutate(
    wave    = factor(wave, levels = WAVES),
    label   = WAVE_LABELS[as.character(wave)],
    lineage = lineage[as.character(wave)],
    cum_incidence_per_100k = round(cases / pop * 1e5, 1),
    cfr_pct = round(deaths / cases * 100, 2)
  ) %>%
  arrange(wave) %>%
  mutate(end = as.Date(dplyr::lead(start) - 1),
         end = as.Date(ifelse(is.na(end), as.Date("2022-05-31"), end),
                       origin = "1970-01-01")) %>%
  select(label, lineage, start, end, n, cases, deaths,
         cum_incidence_per_100k, cfr_pct)

write_csv_audit(wave_tab, "table3_waves.csv")
writexl::write_xlsx(wave_tab, file.path(DIR$tables, "Table3_waves.xlsx"))

# ---- Moran's I -------------------------------------------------------------
# One row per area per wave; neighbours from queen contiguity on the polygons.
one_wave <- dat %>% filter(covid_variant == WAVES[1]) %>% arrange(Etiqueta)
assert(nrow(one_wave) == 82, "expected 82 areas in the first wave")

nb <- poly2nb(one_wave, queen = TRUE)
isolated <- which(card(nb) == 0)
assert(length(isolated) == 0, "some areas have no neighbours under queen contiguity")
listw <- nb2listw(nb, style = "W", zero.policy = TRUE)

moran_one <- function(values) {
  mc <- moran.mc(values, listw, nsim = 1001, zero.policy = TRUE)
  tibble(observed = unname(mc$statistic), p_value = mc$p.value)
}

moran_tab <- lapply(OUTCOMES, function(outcome) {
  var <- switch(outcome,
    incidence = "incidence_rate",
    mortality = "mortality_rate",
    cfr       = "lethality")
  bind_rows(lapply(WAVES, function(w) {
    d <- dat %>% filter(covid_variant == w) %>% arrange(Etiqueta)
    cbind(tibble(outcome = outcome,
                 wave = WAVE_LABELS[[w]]), moran_one(d[[var]]))
  }))
}) %>% bind_rows()

write_csv_audit(moran_tab, "morans_I.csv")
log_msg("Moran's I computed for ", nrow(moran_tab), " outcome-wave combinations")
