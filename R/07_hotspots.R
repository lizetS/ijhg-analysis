# ===========================================================================
# 07_hotspots.R -- hotspot classification, with and without posterior
# uncertainty (Reviewer 2, point 10).
#
# Current rule      : area is elevated in a wave if posterior mean RR > 1.5.
# Uncertainty rule  : area is elevated in a wave if Pr(RR > 1.5 | data) > 0.8
#                     (and > 0.9 as a stricter alternative).
# Persistence is then counted over the seven waves under each rule.
#
# Outputs: output/tables/hotspot_classification.csv
#          output/tables/hotspot_counts.csv
# ===========================================================================


source("config.R")
source(file.path(ROOT, "R", "00_setup.R"))

log_msg("--- 07_hotspots ---")
ar <- utils::read.csv(file.path(DIR$tables, "area_rr_summary.csv"),
                      stringsAsFactors = FALSE) %>%
  filter(outcome == "incidence")

ref <- read_rds("analysis_dataset.rds") %>%
  filter(covid_variant == WAVES[1]) %>%
  st_drop_geometry() %>%
  select(area_code = Etiqueta, ASalud)

count_elevated <- function(df, flag) {
  df %>% group_by(area_code) %>% summarise(n_waves = sum(.data[[flag]]), .groups = "drop")
}

classify <- function(n) {
  case_when(
    n >= CFG$hotspot_min_waves ~ "Persistent",
    n >= 3                     ~ "Frequent",
    n >= 1                     ~ "Occasional",
    TRUE                       ~ "Never"
  )
}

mean_rule  <- ar %>% mutate(elev = rr_mean > CFG$hotspot_rr) %>%
  count_elevated("elev") %>% mutate(rule = "posterior mean RR > 1.5")

prob80_rule <- ar %>% mutate(elev = prob_rr_gt_threshold > 0.80) %>%
  count_elevated("elev") %>% mutate(rule = "Pr(RR > 1.5) > 0.80")

prob90_rule <- ar %>% mutate(elev = prob_rr_gt_threshold > 0.90) %>%
  count_elevated("elev") %>% mutate(rule = "Pr(RR > 1.5) > 0.90")

alls <- bind_rows(mean_rule, prob80_rule, prob90_rule) %>%
  left_join(ref, by = "area_code") %>%
  mutate(category = classify(n_waves))

write_csv_audit(alls, "hotspot_classification.csv")

counts <- alls %>%
  group_by(rule, category) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::complete(rule, category, fill = list(n = 0)) %>%
  mutate(pct = round(100 * n / 82, 1)) %>%
  arrange(rule, factor(category, levels = c("Never", "Occasional", "Frequent", "Persistent")))

write_csv_audit(counts, "hotspot_counts.csv")
print(counts)

# areas reclassified when posterior uncertainty is taken into account
cmp <- alls %>%
  select(area_code, ASalud, rule, category) %>%
  tidyr::pivot_wider(names_from = rule, values_from = category)
reclassified <- cmp %>%
  filter(`posterior mean RR > 1.5` != `Pr(RR > 1.5) > 0.80`)
write_csv_audit(reclassified, "hotspot_reclassified.csv")
log_msg("07_hotspots complete: ", nrow(reclassified),
        " areas change category under the uncertainty rule")
