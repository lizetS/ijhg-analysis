# ===========================================================================
# 23_table2_covariates.R -- Table 2: the covariate set of each model
# (Reviewer 1, point 1.10: list the covariates in each model).
#
# The table is built from covariates_for(), the same function the models use,
# so it cannot drift from the code. It marks, for each of the three outcomes,
# which covariates enter, and flags the two wave-specific ones.
#
# Outputs: output/tables/Table2_covariates.xlsx
#          output/tables/Table2_covariates.csv
# ===========================================================================

source("config.R")
source(file.path(ROOT, "R", "00_setup.R"))

log_msg("--- 23_table2_covariates ---")

# every covariate that can appear in any model, in the order of config.R
all_covs <- unique(c(COV_INCIDENCE, COV_VACCINE, COV_COMORBID, COV_WAVE_VARYING))
outcomes <- c("incidence", "mortality", "cfr")

tab <- tibble(Variable = vapply(all_covs, label_of, character(1)), code = all_covs)
for (oc in outcomes) {
  # the set is the same in every wave except for vaccination coverage, so the
  # Delta model defines the superset
  used <- covariates_for(oc, "Delta")
  tab[[oc]] <- ifelse(all_covs %in% used, "x", "")
}
tab$code <- NULL
tab <- tab %>%
  rename(Incidence = incidence, Mortality = mortality, `Case fatality` = cfr)

notes <- tibble(
  Variable = c(
    "",
    "x = included in the model.",
    "Covariates marked x were rescaled to [0,1] across the 82 health areas; every relative risk is per full observed range.",
    "Vaccination coverage takes a value only in the Delta wave, so it enters the incidence model of that wave alone.",
    "Incidence (same wave) is the crude cumulative incidence of the same wave, log-transformed and standardised within the wave; it is added to the mortality models only, not to case fatality, where cases are already the offset.",
    "Population density is included in the incidence models and omitted from the two severity models: with the six prevalence variables it reaches a variance inflation factor of 5.6 (collinearity.csv)."
  ),
  Incidence = "", Mortality = "", `Case fatality` = ""
)

out <- bind_rows(tab, notes)
write_csv_audit(out, "Table2_covariates.csv")
write_xlsx(out, file.path(DIR$tables, "Table2_covariates.xlsx"))

print(tab)
log_msg("23_table2_covariates complete: ", nrow(tab), " covariables en ",
        length(outcomes), " desenlaces")
