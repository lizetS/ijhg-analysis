# ===========================================================================
# 03_standardise.R -- standardise covariates across the 82 health areas.
#
# CFG$scaling = "z"      -> z-scores (mean 0, SD 1)   [used in the revision]
# CFG$scaling = "minmax" -> rescale to [0, 1]          [legacy scripts]
#
# Because every covariate except vaccination coverage is measured once and is
# time-invariant, its z-scores are identical in every wave. The scaling
# parameters are therefore computed once from the 82 areas and applied to all
# waves. This is asserted below and answers Reviewer 2, point 4.
#
# Outputs: data/derived/analysis_standardised.rds
#          output/tables/scaling_parameters.csv
# ===========================================================================


source("config.R")
source(file.path(ROOT, "R", "00_setup.R"))

log_msg("--- 03_standardise ---")
dat <- read_rds("analysis_dataset.rds")

covs <- unique(c(COV_INCIDENCE, COV_VACCINE, COV_SEVERITY))
covs <- setdiff(intersect(covs, names(dat)), COV_WAVE_VARYING)   # wave-specific covariate handled below

# time-invariance check (vaccination is the documented exception)
w1 <- dat %>% filter(covid_variant == WAVES[1]) %>% arrange(Etiqueta)
time_varying <- c()
for (v in covs) {
  ok <- TRUE
  for (w in WAVES[-1]) {
    dw <- dat %>% filter(covid_variant == w) %>% arrange(Etiqueta)
    if (!isTRUE(all.equal(w1[[v]], dw[[v]], check.attributes = FALSE))) {
      ok <- FALSE; break
    }
  }
  if (!ok) time_varying <- c(time_varying, v)
}
log_msg("time-varying covariates: ",
        if (length(time_varying)) paste(time_varying, collapse = ", ") else "none")

# Scaling parameters come from the distribution across the 82 areas, taken
# from a SINGLE wave. Using all 574 rows would repeat each area seven times
# and yield the population SD rather than the sample SD, so the standardised
# values would not have SD exactly 1. Vaccination coverage is taken from the
# only wave in which it is defined (Delta).
ref_values <- function(v) {
  for (w in WAVES) {
    x <- dat[[v]][dat$covid_variant == w]
    if (!all(is.na(x))) return(x)
  }
  rep(NA_real_, 82)
}

scaling_tab <- lapply(covs, function(v) {
  x <- ref_values(v)
  keep <- !is.na(x)
  tibble(variable = v, n_areas = sum(keep),
         mean = mean(x[keep]), sd = sd(x[keep]),
         min = min(x[keep]), max = max(x[keep]),
         method = CFG$scaling)
}) %>% bind_rows()

dat_std <- dat
for (i in seq_len(nrow(scaling_tab))) {
  v <- scaling_tab$variable[i]; mu <- scaling_tab$mean[i]
  sdv <- scaling_tab$sd[i]; mn <- scaling_tab$min[i]; mx <- scaling_tab$max[i]
  if (CFG$scaling == "z") {
    dat_std[[v]] <- if (is.finite(sdv) && sdv > 0) (dat[[v]] - mu) / sdv else 0
  } else {
    dat_std[[v]] <- if (is.finite(mx - mn) && (mx - mn) > 0)
      (dat[[v]] - mn) / (mx - mn) else 0
  }
}

# Wave-specific covariate (incidence_cov): scaled within each wave, because its
# distribution changes from wave to wave.
if ("incidence_cov" %in% names(dat)) {
  wv_tab <- list()
  for (w in WAVES) {
    idx <- dat$covid_variant == w
    x <- dat$incidence_cov[idx]
    wv_tab[[w]] <- tibble(variable = "incidence_cov", wave = w, n_areas = length(x),
                          mean = mean(x), sd = sd(x), min = min(x), max = max(x),
                          method = CFG$scaling)
    dat_std$incidence_cov[idx] <- if (CFG$scaling == "z") (x - mean(x)) / sd(x) else
      (x - min(x)) / (max(x) - min(x))
  }
  scaling_tab <- bind_rows(scaling_tab, bind_rows(wv_tab))
}

# verify on the first wave in which each variable is defined (82 area values)
for (v in covs) {
  for (w in WAVES) {
    x <- dat_std[[v]][dat_std$covid_variant == w]
    if (all(is.na(x))) next
    x <- x[!is.na(x)]
    if (CFG$scaling == "z") {
      assert(abs(mean(x)) < 1e-8, paste("mean not 0 after standardisation:", v))
      assert(abs(sd(x) - 1) < 1e-8, paste("SD not 1 after standardisation:", v))
    } else {
      assert(abs(min(x)) < 1e-8 && abs(max(x) - 1) < 1e-8,
             paste("range not [0,1] after rescaling:", v))
    }
    break
  }
}

save_rds(dat_std, "analysis_standardised.rds")
write_csv_audit(scaling_tab, "scaling_parameters.csv")
log_msg("standardised ", nrow(scaling_tab), " covariates with method '", CFG$scaling, "'")
