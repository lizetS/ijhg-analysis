# ===========================================================================
# 06_null_interval.R -- Supplementary Table S3.
#
# Reclassifies every estimate under three criteria:
#   primary : 95% CrI excludes 1.0            (reported criterion)
#   band    : 95% CrI excludes 0.9-1.1        (original, asymmetric)
#   sym     : 95% CrI excludes 0.8-1.25       (symmetric, stricter)
#
# Output: output/tables/TableS3_null_interval.csv
# ===========================================================================


source("config.R")
source(file.path(ROOT, "R", "00_setup.R"))

log_msg("--- 06_null_interval ---")
co <- utils::read.csv(file.path(DIR$tables, "coefficients.csv"),
                      stringsAsFactors = FALSE) %>%
  filter(variable != "(Intercept)")

excludes <- function(lo, hi, band) hi < band[1] | lo > band[2]

s3 <- co %>%
  mutate(
    excludes_primary = excludes(lower, upper, CFG$null_primary),
    excludes_band    = excludes(lower, upper, CFG$null_band),
    excludes_sym     = excludes(lower, upper, CFG$null_sym)
  ) %>%
  select(outcome, wave, variable, rr, lower, upper,
         excludes_primary, excludes_band, excludes_sym) %>%
  mutate(variable = vapply(variable, label_of, character(1)))

write_csv_audit(s3, "TableS3_null_interval.csv")

summary_tab <- s3 %>%
  summarise(
    n_estimates      = n(),
    n_primary        = sum(excludes_primary),
    n_band           = sum(excludes_band),
    n_sym            = sum(excludes_sym),
    pct_primary      = round(100 * mean(excludes_primary), 1),
    pct_band         = round(100 * mean(excludes_band), 1),
    pct_sym          = round(100 * mean(excludes_sym), 1)
  )
write_csv_audit(summary_tab, "TableS3_summary.csv")
print(summary_tab)

# Estimates that qualify under the band but not under the stricter symmetric
# criterion - the ones whose status is criterion-dependent.
discordant <- s3 %>% filter(excludes_band & !excludes_sym)
write_csv_audit(discordant, "TableS3_criterion_dependent.csv")
log_msg("06_null_interval complete: ", nrow(s3), " estimates, ",
        nrow(discordant), " depend on the criterion")
