# ===========================================================================
# 04_fit_models.R -- fit the wave-specific Leroux CAR models.
#
# One model per outcome x wave. THREE chains per model via CARBayes' native
# n.chains argument, combined into an mcmc.list, with Gelman-Rubin,
# effective sample size and ACF diagnostics (Reviewer 2, point 6).
# Area-level posterior summaries are computed inside the fit so that hotspot
# classification (point 10) can use posterior probabilities.
#
# Outputs
#   output/tables/coefficients.csv        tidy RR table for all models
#   output/tables/diagnostics.csv         Rhat / ESS / status per model
#   output/tables/area_rr_summary.csv     posterior mean and Pr(RR > 1.5)
#   output/models/<outcome>_<wave>.rds    compact model object
# ===========================================================================


source("config.R")
source(file.path(ROOT, "R", "00_setup.R"))

log_msg("--- 04_fit_models ---")
DAT <- read_rds("analysis_standardised.rds")

DF <- DAT %>% st_drop_geometry() %>% arrange(Etiqueta, covid_variant)

SP <- DAT %>% filter(covid_variant == WAVES[1]) %>% arrange(Etiqueta)
NB <- poly2nb(SP, queen = TRUE)
W_QUEEN <- nb2mat(NB, style = "B")
assert(nrow(W_QUEEN) == 82 && all(rowSums(W_QUEEN) > 0),
       "adjacency matrix has isolated areas")

# ---------------------------------------------------------------------------
# Core fitter. Used by this script and by 08_sensitivity.R.
#   outcome : "incidence" | "mortality" | "cfr"
#   wave    : element of WAVES
#   covs    : character vector of covariates
#   W       : adjacency matrix
#   family  : "poisson" (default) or "binomial" (CFR sensitivity)
# ---------------------------------------------------------------------------
fit_model <- function(outcome, wave, covs, W = W_QUEEN,
                      family = "poisson", priors = list()) {
  d <- DF %>% filter(covid_variant == wave) %>% arrange(Etiqueta)

  dropped <- 0L
  if (outcome == "cfr" || family == "binomial") {
    keep <- d$cases_per_variant > 0
    dropped <- sum(!keep)
    d <- d[keep, , drop = FALSE]
  }

  f <- if (family == "binomial") {
    stats::as.formula(paste("deaths_per_variant ~",
                            paste(covs, collapse = " + ")))
  } else {
    build_formula(covs, outcome)
  }

  X <- model.matrix(f, data = d)
  areas <- d$Etiqueta

  if (outcome == "incidence") {
    if (identical(Sys.getenv("QUICK_RUN"), "1") || nzchar(Sys.getenv("ITER_OVERRIDE"))) {
      n_it <- CFG$n_iter_incidence; b_in <- CFG$burnin_incidence; th <- CFG$thin_incidence
    } else {
      row <- INC_PLAN[INC_PLAN$wave == wave, ]
      n_it <- as.integer(row$n_iter); b_in <- as.integer(row$burnin); th <- as.integer(row$thin)
    }
  } else {
    n_it <- CFG$n_iter_severity;  b_in <- CFG$burnin_severity;  th <- CFG$thin_severity
  }

  set.seed(CFG$seed + match(wave, WAVES))
  args <- list(formula = f, family = family, data = d, W = W,
               burnin = b_in, n.sample = n_it, thin = th,
               n.chains = CFG$n_chains, n.cores = 1, verbose = FALSE)
  if (family == "binomial") args$trials <- d$cases_per_variant
  if (!is.null(priors$prior.mean.beta)) args$prior.mean.beta <- priors$prior.mean.beta
  if (!is.null(priors$prior.var.beta))  args$prior.var.beta  <- priors$prior.var.beta
  if (!is.null(priors$prior.tau2))      args$prior.tau2      <- priors$prior.tau2

  # Very sparse outcomes (Wuhan mortality: 20 deaths, 66 areas at zero) make the
  # MALA acceptance ratio NaN and abort CARBayes. Fall back to the random-walk
  # sampler, then to a tighter beta prior, and record what was used.
  run_sampler <- function(mala, pvb) {
    args$MALA <- mala
    if (!is.null(pvb)) args$prior.var.beta <- pvb
    do.call(S.CARleroux, args)
  }
  tag <- paste(outcome, wave, sep = "_")
  sampler <- "MALA"
  args$MALA <- TRUE
  mod <- tryCatch(run_sampler(TRUE, NULL), error = function(e1) {
    log_msg("MALA failed for ", tag, " (", conditionMessage(e1),
            "); retrying with MALA = FALSE")
    sampler <<- "random walk"
    tryCatch(run_sampler(FALSE, NULL), error = function(e2) {
      if (!is.null(priors$prior.var.beta)) stop(e2)
      log_msg("random walk failed too; retrying with prior.var.beta = 10")
      sampler <<- "random walk, prior.var.beta = 10"
      run_sampler(FALSE, rep(10, ncol(X)))
    })
  })

  # CARBayes returns a named list of parameter blocks, each an mcmc.list with
  # one element per chain.
  beta_ml <- mod$samples$beta
  phi_ml  <- mod$samples$phi
  assert(length(beta_ml) == CFG$n_chains,
         "CARBayes returned an unexpected number of chains")

  gd  <- gelman.diag(beta_ml, autoburnin = FALSE, multivariate = FALSE)
  ess <- effectiveSize(beta_ml)
  pooled <- as.matrix(beta_ml)

  coef_names <- c("(Intercept)", covs)
  assert(ncol(pooled) == length(coef_names), "unexpected number of beta parameters")
  colnames(pooled) <- coef_names

  coef_tab <- tibble(
    outcome = outcome, wave = WAVE_LABELS[[wave]],
    variable = coef_names,
    rr = exp(colMeans(pooled)),
    lower = exp(apply(pooled, 2, quantile, 0.025)),
    upper = exp(apply(pooled, 2, quantile, 0.975)),
    pr_gt_1 = colMeans(pooled > 0)
  )

  # area-level posterior summaries, accumulated in chunks
  # exp(X beta + phi) is a rate per person, not a relative risk. The area RR
  # is expressed relative to the population-weighted mean rate of the city,
  # so that RR = 1 means "at the city average" and thresholds such as 1.5 are
  # interpretable.
  n <- nrow(X)
  w_pop <- d$poblacion
  acc_mean <- numeric(n); acc_prob <- numeric(n); total <- 0
  for (ch in seq_along(phi_ml)) {
    B   <- as.matrix(beta_ml[[ch]])
    Phi <- as.matrix(phi_ml[[ch]])
    assert(ncol(B) == ncol(X), "beta/design mismatch")
    assert(ncol(Phi) == n, "phi/area mismatch")
    idx <- seq_len(nrow(B))
    for (k in split(idx, ceiling(idx / 5000))) {
      lam <- exp(B[k, , drop = FALSE] %*% t(X) + Phi[k, , drop = FALSE])
      overall <- as.numeric(lam %*% w_pop) / sum(w_pop)
      RR <- lam / overall
      acc_mean <- acc_mean + colSums(RR)
      acc_prob <- acc_prob + colSums(RR > CFG$hotspot_rr)
      total <- total + length(k)
    }
  }
  area_tab <- tibble(area_code = areas,
                     rr_mean = acc_mean / total,
                     prob_rr_gt_threshold = acc_prob / total)

  rho_mean  <- mean(vapply(mod$samples$rho,  function(S) mean(as.matrix(S)[, 1]), numeric(1)))
  tau2_mean <- mean(vapply(mod$samples$tau2, function(S) mean(as.matrix(S)[, 1]), numeric(1)))

  rhat_max <- max(gd$psrf[, "Point est."], na.rm = TRUE)
  ess_min  <- min(ess)
  conv_status <- if (is.finite(rhat_max) && rhat_max < 1.05 &&
                     is.finite(ess_min) && ess_min >= 100) "ok" else "check convergence"

  # Pearson dispersion (Reviewer 2, point 7). Plug-in version based on the
  # posterior-mean fitted values: X2 = sum((y - mu)^2 / var(mu)), dispersion =
  # X2 / (n - p) with p the number of regression coefficients (the random
  # effects also absorb degrees of freedom, so X2 / n is reported too).
  # ppp = share of replicated data sets (drawn from the same likelihood at the
  # fitted means) with a larger X2 than the observed one; values near 0 mean
  # more variability than the likelihood allows (overdispersion).
  disp <- tryCatch({
    mu <- as.numeric(mod$fitted.values)
    y  <- if (outcome == "incidence") d$cases_per_variant else d$deaths_per_variant
    ok <- is.finite(mu) & mu > 1e-8 & length(mu) == length(y)
    assert(all(ok), "fitted values not usable for the Pearson statistic")
    if (family == "binomial") {
      nt  <- d$cases_per_variant
      v   <- mu * (1 - mu / nt)
      rep_y <- function() stats::rbinom(length(mu), nt, mu / nt)
    } else {
      v   <- mu
      rep_y <- function() stats::rpois(length(mu), mu)
    }
    keep <- v > 1e-8
    x2   <- sum((y[keep] - mu[keep])^2 / v[keep])
    x2r  <- vapply(seq_len(1000), function(i) {
      yr <- rep_y(); sum((yr[keep] - mu[keep])^2 / v[keep])
    }, numeric(1))

    # Moran's I of the Pearson residuals (Reviewer 2, point 6): does the CAR
    # term leave spatial structure in what the model cannot explain? The
    # residuals follow the area order of W, which is the order of d (both are
    # arranged by Etiqueta).
    moran <- c(stat = NA_real_, p = NA_real_)
    if (length(y) == nrow(W) && all(keep)) {
      res <- (y - mu) / sqrt(v)
      lw  <- spdep::mat2listw(W, style = "W")
      mt  <- spdep::moran.test(res, lw, zero.policy = TRUE)
      moran <- c(stat = unname(mt$estimate["Moran I statistic"]), p = mt$p.value)
    }

    list(chi2 = x2, disp = x2 / (length(y) - length(coef_names)),
         per_obs = x2 / length(y), ppp = mean(x2r >= x2),
         moran = moran["stat"], moran_p = moran["p"])
  }, error = function(e) {
    log_msg("  Pearson statistic not computed: ", conditionMessage(e))
    list(chi2 = NA_real_, disp = NA_real_, per_obs = NA_real_, ppp = NA_real_,
         moran = NA_real_, moran_p = NA_real_)
  })

  # Computed before the tibble: inside it the earlier column
  # `wave = WAVE_LABELS[[wave]]` shadows the argument, so comparing
  # DF$covid_variant to `wave` there matched nothing for Wuhan1/Wuhan2 (whose
  # labels are "Early"/"Wuhan") and n_deaths came out as 0.
  nd_wave <- sum(DF$deaths_per_variant[DF$covid_variant == wave])

  diag_tab <- tibble(
    outcome = outcome, wave = WAVE_LABELS[[wave]],
    n_areas = nrow(d), n_dropped_zero_offset = dropped,
    n_chains = CFG$n_chains, n_iter = n_it,
    burnin = b_in, thin = th, sampler = sampler,
    rhat_max = rhat_max, ess_min = ess_min,
    rho_mean = rho_mean, tau2_mean = tau2_mean,
    pearson_chi2 = disp$chi2, pearson_dispersion = disp$disp,
    pearson_chi2_per_obs = disp$per_obs, pearson_ppp = disp$ppp,
    moran_residual = unname(disp$moran), moran_p = unname(disp$moran_p),
    n_deaths = nd_wave,
    sparse_outcome = outcome != "incidence" && nd_wave < MIN_DEATHS_SEVERITY,
    status = conv_status
  )

  keep_idx <- seq(1, nrow(pooled), by = 100)
  chain_id <- rep(seq_len(CFG$n_chains), each = nrow(pooled) / CFG$n_chains)
  obj <- list(
    outcome = outcome, wave = wave, formula = f, covariates = covs,
    coefficients = coef_tab, diagnostics = diag_tab, areas = area_tab,
    trace = cbind(as.data.frame(pooled[keep_idx, , drop = FALSE]),
                  chain = chain_id[keep_idx]),
    rho_trace = unlist(lapply(mod$samples$rho, function(S)
      as.matrix(S)[seq(1, nrow(as.matrix(S)), by = 100), 1])),
    config = CFG
  )

  list(obj = obj, coef = coef_tab, diag = diag_tab, areas = area_tab)
}

# ---------------------------------------------------------------------------
# Run the primary models
# ---------------------------------------------------------------------------
if (!exists("RUN_SENSITIVITY")) {
  results <- list()
  for (outcome in OUTCOMES) {
    for (wave in WAVES) {
      tag <- paste(outcome, wave, sep = "_")

      # severity models: first five waves only (Omicron left out, see config.R)
      if (outcome != "incidence") {
        nd <- sum(DF$deaths_per_variant[DF$covid_variant == wave])
        if (!(wave %in% SEVERITY_WAVES)) {
          log_msg("skipping ", tag, ": wave not in SEVERITY_WAVES (", nd, " deaths)")
          results[[tag]] <- list(
            obj = NULL, coef = NULL, areas = NULL,
            diag = tibble(outcome = outcome, wave = WAVE_LABELS[[wave]],
                          n_deaths = nd,
                          status = sprintf("not estimated (excluded wave, %d deaths)", nd)))
          next
        }
        if (nd < MIN_DEATHS_SEVERITY) {
          log_msg("WARNING ", tag, ": sparse outcome, only ", nd, " deaths")
        }
      }

      log_msg("fitting ", tag)
      res <- tryCatch(fit_model(outcome, wave, covariates_for(outcome, wave)),
                      error = function(e) {
                        log_msg("FAILED ", tag, ": ", conditionMessage(e))
                        list(obj = NULL, coef = NULL,
                             diag = tibble(outcome = outcome,
                                           wave = WAVE_LABELS[[wave]],
                                           status = paste("failed:", conditionMessage(e))),
                             areas = NULL)
                      })
      results[[tag]] <- res
      if (!is.null(res$obj)) save_rds(res$obj, paste0(tag, ".rds"), DIR$models)
    }
  }

  coef_all <- bind_rows(lapply(results, `[[`, "coef"))
  diag_all <- bind_rows(lapply(results, `[[`, "diag"))
  area_all <- bind_rows(lapply(results, function(r) {
    if (is.null(r$areas)) return(NULL)
    cbind(tibble(outcome = r$obj$outcome, wave = WAVE_LABELS[[r$obj$wave]]), r$areas)
  }))

  write_csv_audit(coef_all, "coefficients.csv")
  write_csv_audit(diag_all, "diagnostics.csv")
  write_csv_audit(area_all, "area_rr_summary.csv")

  n_ok  <- sum(diag_all$status == "ok")
  n_net <- sum(grepl("^not estimated", diag_all$status))
  n_fail <- nrow(diag_all) - n_ok - n_net
  log_msg("04_fit_models complete: ", n_ok, " fitted, ",
          n_net, " not estimated (event count), ", n_fail, " failed")
  print(diag_all %>% select(outcome, wave, rhat_max, ess_min, status))
}
