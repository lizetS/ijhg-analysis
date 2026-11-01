# ===========================================================================
# 05_tables.R -- manuscript Table 1 and Supplementary Tables S1/S2.
#
# Rules applied (all configurable in config.R):
#   bold-style flag : 95% CrI excludes 1.0            (CFG$null_primary)
#   dagger          : interval ratio upper/lower > 20  (sparse / unstable)
#
# Outputs: output/tables/Table1_incidence.xlsx
#          output/tables/TableS1_mortality.xlsx
#          output/tables/TableS2_case_fatality.xlsx
#          output/tables/table1_long.csv (machine-readable)
# ===========================================================================


source("config.R")
source(file.path(ROOT, "R", "00_setup.R"))

log_msg("--- 05_tables ---")
co <- utils::read.csv(file.path(DIR$tables, "coefficients.csv"),
                      stringsAsFactors = FALSE)

co <- co %>% filter(variable != "(Intercept)")

WIDE_RATIO <- 20

co <- co %>%
  mutate(
    excludes_null = lower > CFG$null_primary[1] | upper < CFG$null_primary[2],
    sparse        = (upper / lower) > WIDE_RATIO,
    cell          = sprintf("%.2f (%.2f–%.2f)%s",
                            rr, lower, upper, ifelse(sparse, " †", ""))
  )

build_wide <- function(outcome) {
  d <- co %>% filter(outcome == !!outcome)
  d$variable <- factor(d$variable, levels = unique(d$variable))
  wide <- d %>%
    select(variable, wave, cell) %>%
    tidyr::pivot_wider(names_from = wave, values_from = cell)
  wide
}

tables <- list(
  incidence = build_wide("incidence"),
  mortality = build_wide("mortality"),
  cfr       = build_wide("cfr")
)

for (nm in names(tables)) {
  tab <- tables[[nm]]
  tab$Variable <- vapply(as.character(tab$variable), label_of, character(1))
  tab <- tab %>% select(Variable, everything(), -variable)
  file <- switch(nm,
    incidence = "Table1_incidence.xlsx",
    mortality = "TableS1_mortality.xlsx",
    cfr       = "TableS2_case_fatality.xlsx")
  path <- file.path(DIR$tables, file)
  write_xlsx(tab, path)
  log_msg("wrote ", path)
}

write_csv_audit(co %>% select(outcome, wave, variable, rr, lower, upper,
                              pr_gt_1, excludes_null, sparse),
                "table1_long.csv")

log_msg("05_tables complete")
