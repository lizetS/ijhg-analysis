# ===========================================================================
# model_worker.R -- fit a subset of the models in its own process.
#
# Launched by tools/run_parallel.R, one process per worker. Each model is
# seeded inside fit_model from (seed + wave), so a model gets the same draws
# whatever worker fits it: splitting the work does not change the results.
#
#   PROJ    repository folder
#   MODELS  comma-separated outcome:wave pairs, e.g. "incidence:Delta,cfr:Beta"
#   OUT     folder for the partial coefficient/diagnostic/area files
# ===========================================================================

# Parameters come as command-line arguments: on Windows system2(wait = FALSE)
# does not propagate the `env` argument, so environment variables are not
# reliable here. Fall back to them only if no arguments were given.
a <- commandArgs(trailingOnly = TRUE)
proj <- if (length(a) >= 1) a[1] else Sys.getenv("PROJ")
mods <- if (length(a) >= 2) a[2] else Sys.getenv("MODELS")
out  <- if (length(a) >= 3) a[3] else Sys.getenv("OUT")
stopifnot(nzchar(proj), nzchar(mods), nzchar(out))

# A variant run (alternative crowding source, replication spec) is selected by
# RUN_TAG / HACIN_SOURCE, and config.R reads them from the environment. They are
# passed as arguments because system2(wait = FALSE) does not forward `env` on
# Windows, so they have to be set here, before config.R is sourced.
if (length(a) >= 4 && nzchar(a[4])) Sys.setenv(RUN_TAG = a[4])
if (length(a) >= 5 && nzchar(a[5])) Sys.setenv(HACIN_SOURCE = a[5])
setwd(proj)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

source("config.R")
source(file.path(ROOT, "R", "00_setup.R"))
RUN_SENSITIVITY <- TRUE
source(file.path(ROOT, "R", "04_fit_models.R"))   # defines fit_model, DF, W_QUEEN

for (p in strsplit(mods, ",", fixed = TRUE)[[1]]) {
  parts <- strsplit(p, ":", fixed = TRUE)[[1]]
  oc <- parts[1]; wv <- parts[2]
  tag <- paste(oc, wv, sep = "_")
  if (file.exists(file.path(out, paste0(tag, ".done")))) next

  log_msg("[worker ", Sys.getpid(), "] fitting ", tag)
  t0 <- Sys.time()
  r <- tryCatch(fit_model(oc, wv, covariates_for(oc, wv)),
                error = function(e) {
                  log_msg("[worker ", Sys.getpid(), "] FAILED ", tag, ": ",
                          conditionMessage(e))
                  NULL
                })

  if (is.null(r)) {
    utils::write.csv(
      tibble(outcome = oc, wave = WAVE_LABELS[[wv]], status = "failed"),
      file.path(out, paste0(tag, "_diag.csv")), row.names = FALSE)
  } else {
    saveRDS(r$obj, file.path(DIR$models, paste0(tag, ".rds")))
    utils::write.csv(r$coef, file.path(out, paste0(tag, "_coef.csv")), row.names = FALSE)
    utils::write.csv(r$diag, file.path(out, paste0(tag, "_diag.csv")), row.names = FALSE)
    utils::write.csv(
      cbind(tibble(outcome = oc, wave = WAVE_LABELS[[wv]]), r$areas),
      file.path(out, paste0(tag, "_area.csv")), row.names = FALSE)
    log_msg("[worker ", Sys.getpid(), "] done ", tag, " in ",
            round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
  }
  file.create(file.path(out, paste0(tag, ".done")))
  rm(r); invisible(gc())
}

log_msg("[worker ", Sys.getpid(), "] finished")
