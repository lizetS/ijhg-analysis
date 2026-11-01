# ===========================================================================
# 26_bym.R -- Poisson models with an additional unstructured random effect
# (BYM) for incidence and mortality, as promised in response 2.7.
#
# S.CARbym = spatial Leroux-type component + independent area effect. It does
# not have a negative-binomial likelihood; this complements 24_negative_binomial
# (non-spatial NB) by absorbing extra-Poisson variation at area level.
#
# Output: output/tables/sensitivity_bym.csv
# ===========================================================================

RUN_SENSITIVITY <- TRUE
source("R/04_fit_models.R")     # DF, SP, W_QUEEN, CFG, build_formula

log_msg("--- 26_bym ---")

fit_bym <- function(oc, wave) {
  d <- DF %>% filter(covid_variant == wave) %>% arrange(Etiqueta)
  covs <- intersect(covariates_for(oc, wave), names(d))
  f <- build_formula(covs, oc)
  if (oc == "incidence") {
    row <- INC_PLAN[INC_PLAN$wave == wave, ]
    n_it <- as.integer(row$n_iter); b_in <- as.integer(row$burnin); th <- as.integer(row$thin)
  } else {
    n_it <- CFG$n_iter_severity; b_in <- CFG$burnin_severity; th <- CFG$thin_severity
  }
  set.seed(CFG$seed + match(wave, WAVES) + 10000L)
  mod <- tryCatch(CARBayes::S.CARbym(
    formula = f, family = "poisson", data = d, W = W_QUEEN,
    burnin = b_in, n.sample = n_it, thin = th, n.chains = CFG$n_chains,
    n.cores = 1, verbose = FALSE), error = function(e) {
      log_msg("  BYM failed for ", oc, "_", wave, ": ", conditionMessage(e)); NULL
    })
  if (is.null(mod)) return(NULL)
  bm <- mod$samples$beta
  pooled <- as.matrix(bm); colnames(pooled) <- c("(Intercept)", covs)
  tibble(outcome = oc, wave = WAVE_LABELS[[wave]], variable = colnames(pooled),
         rr = exp(colMeans(pooled)),
         lower = exp(apply(pooled, 2, stats::quantile, 0.025)),
         upper = exp(apply(pooled, 2, stats::quantile, 0.975)),
         pr_gt_1 = colMeans(pooled > 0),
         model = "BYM")
}

jobs <- rbind(
  data.frame(oc = "incidence", wave = WAVES, stringsAsFactors = FALSE),
  data.frame(oc = "mortality", wave = SEVERITY_WAVES, stringsAsFactors = FALSE))

res <- parallel::mclapply(seq_len(nrow(jobs)), function(i)
  fit_bym(jobs$oc[i], jobs$wave[i]), mc.cores = 6)

out <- bind_rows(res) %>% filter(variable != "(Intercept)")
write_csv_audit(out, "sensitivity_bym.csv")
log_msg("26_bym complete: ", nrow(out), " rows")
