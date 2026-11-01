# ===========================================================================
# sens_worker.R -- fit a subset of the sensitivity models in its own process.
#
# Launched by tools/run_sens_parallel.R. Each job is a spec
#   <section>;<outcome>;<wave>;<variant>
# where section is A (comorbidities), A2 (reduced set), B (binomial CFR),
# C (alternative spatial weights) or D (alternative CAR priors). The variant
# names the weight matrix (queen, rook, queen2, knn4) or the prior set
# (default, ig_1_001, ig_half); "-" means not applicable.
#
# Writes one coefficients file per job into OUT and a .done marker, so the
# driver can resume and merge without holding any model in memory.
# ===========================================================================

a <- commandArgs(trailingOnly = TRUE)
proj <- a[1]; jobspec <- a[2]; out <- a[3]
stopifnot(nzchar(proj), nzchar(jobspec), nzchar(out))
setwd(proj)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

source("config.R")
source(file.path(ROOT, "R", "00_setup.R"))
RUN_SENSITIVITY <- TRUE
source(file.path(ROOT, "R", "04_fit_models.R"))   # fit_model, DF, SP, W_QUEEN, NB

# --- alternative spatial weights, built once per worker ---------------------
coords <- st_coordinates(st_centroid(SP))
W_alt <- list(
  queen  = W_QUEEN,
  rook   = nb2mat(poly2nb(SP, queen = FALSE), style = "B"),
  queen2 = nb2mat(nblag_cumul(nblag(NB, 2L)), style = "B", zero.policy = TRUE),
  knn4   = nb2mat(knn2nb(knearneigh(coords, k = 4)), style = "B")
)
# CARBayes requires a symmetric W; the k-NN graph is symmetrised.
W_alt <- lapply(W_alt, function(M) { M <- 1L * ((M + t(M)) > 0); M })  # symmetrise k-NN

# CARBayes' default for tau2 IS c(1, 0.01), so a set of c(1, 0.01) would be the
# same prior twice. The three sets here are genuinely distinct.
prior_sets <- list(
  default = list(),
  ig_half = list(prior.tau2 = c(0.5, 0.0005)),
  ig_vague = list(prior.tau2 = c(0.001, 0.001))
)

for (spec in strsplit(jobspec, ",", fixed = TRUE)[[1]]) {
  f <- strsplit(spec, ";", fixed = TRUE)[[1]]
  sec <- f[1]; oc <- f[2]; wv <- f[3]; var <- if (length(f) >= 4) f[4] else "-"
  tag <- paste(sec, oc, wv, var, sep = "_")
  if (file.exists(file.path(out, paste0(tag, ".done")))) next

  covs <- switch(sec,
    A  = c(covariates_for("incidence", wv), COV_COMORBID),
    A2 = c(COV_INCIDENCE_REDUCED, if (wv == "Delta") COV_VACCINE),
    B  = covariates_for("cfr", wv),
    C  = covariates_for("incidence", wv),
    D  = covariates_for("incidence", wv))
  covs <- intersect(covs, names(DF))

  W <- if (sec == "C") W_alt[[var]] else W_QUEEN
  pr <- if (sec == "D") prior_sets[[var]] else list()
  fam <- if (sec == "B") "binomial" else "poisson"

  log_msg("[sens ", Sys.getpid(), "] ", tag)
  t0 <- Sys.time()
  r <- tryCatch(fit_model(oc, wv, covs, W = W, family = fam, priors = pr),
                error = function(e) { log_msg("  FAILED ", tag, ": ",
                                              conditionMessage(e)); NULL })
  if (!is.null(r)) {
    utils::write.csv(r$coef, file.path(out, paste0(tag, ".csv")), row.names = FALSE)
    log_msg("[sens ", Sys.getpid(), "] done ", tag, " in ",
            round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
  }
  file.create(file.path(out, paste0(tag, ".done")))
  rm(r); invisible(gc())
}
log_msg("[sens ", Sys.getpid(), "] finished")
