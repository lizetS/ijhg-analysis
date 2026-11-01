# ===========================================================================
# run_parallel.R -- fit every model in N worker processes, then merge.
#
# Usage (from the repository folder):
#   N_WORKERS=3 Rscript tools/run_parallel.R
#
# Workers are balanced by estimated cost (iterations): the longest chain goes
# first, and each model goes to the least-loaded worker. Each model is seeded
# from (seed + wave) inside fit_model, so results do not depend on scheduling.
# Partial files go to output/tables/_workers/ and are merged into
# coefficients.csv, diagnostics.csv and area_rr_summary.csv at the end.
#
# Downstream scripts are not run here; run 05, 06, 07, 09, 16, 19, 17 and 08c
# afterwards, or call run_all.R for the rest.
# ===========================================================================

proj <- Sys.getenv("PROJ", unset = getwd())
a    <- commandArgs(trailingOnly = TRUE)
N    <- if (length(a) >= 1) as.integer(a[1]) else
  as.integer(Sys.getenv("N_WORKERS", "2"))
RESUME <- identical(if (length(a) >= 2) a[2] else Sys.getenv("RESUME", "0"), "1")
# optional variant selectors, forwarded to the workers as arguments
RUN_TAG_ARG    <- if (length(a) >= 3) a[3] else ""
HACIN_SRC_ARG  <- if (length(a) >= 4) a[4] else ""
if (nzchar(RUN_TAG_ARG))   Sys.setenv(RUN_TAG = RUN_TAG_ARG)
if (nzchar(HACIN_SRC_ARG)) Sys.setenv(HACIN_SOURCE = HACIN_SRC_ARG)
stopifnot(!is.na(N), N >= 1)

setwd(proj)
source("config.R")
source(file.path(ROOT, "R", "00_setup.R"))

# Refuse to start while another run is writing models in this project: two runs
# writing output/models and output/tables at once corrupt each other's results
# and exhaust memory. PLAN_ONLY only writes the plan, so it is always allowed.
# FORCE_RUN=1 skips the "another run is active" guard, for a deliberate refit
# started right after a previous pass wrote its models.
if (!identical(Sys.getenv("PLAN_ONLY"), "1") &&
    !identical(Sys.getenv("FORCE_RUN"), "1")) {
  # A resume run is expected to start right after the previous one wrote models,
  # so only the log check applies then.
  recent <- list.files(DIR$models, pattern = "\\.rds$", full.names = TRUE)
  fresh_model <- !RESUME && length(recent) &&
    any(as.numeric(difftime(Sys.time(), file.mtime(recent), units = "mins")) < 15)

  # A run is also active when the last "fitting" line has no matching model
  # written after it: the model is still in progress. The age check keeps a
  # crashed run from blocking the launcher for ever.
  logf <- file.path(DIR$logs, "run.log")
  fitting_now <- FALSE
  if (file.exists(logf)) {
    lines <- readLines(logf, warn = FALSE)
    # only the driver's own workers count: scripts 08 and 20 also log "fitting"
    # but write no model, and their dangling lines used to block the launcher
    fit_i   <- suppressWarnings(max(grep("\\[worker .*\\] fitting ", lines), na.rm = TRUE))
    wrote_i <- suppressWarnings(max(grep("models.*\\.rds", lines), na.rm = TRUE))
    # a run that logged its final line is over, whatever the last model line says
    done_i  <- suppressWarnings(max(grep("run_parallel complete", lines), na.rm = TRUE))
    fresh_log <- as.numeric(difftime(Sys.time(), file.mtime(logf), units = "mins")) < 180
    if (is.finite(fit_i) && (!is.finite(wrote_i) || fit_i > wrote_i) &&
        (!is.finite(done_i) || done_i < fit_i) && fresh_log)
      fitting_now <- TRUE
  }

  if (fresh_model || fitting_now) {
    log_msg("run_parallel: another run is active in this project. Refusing to ",
            "start; stop the other run first.")
    quit(save = "no", status = 2)
  }
}

# ---- the model list and its estimated cost --------------------------------
jobs <- list()
for (oc in OUTCOMES) {
  for (w in WAVES) {
    if (oc != "incidence" && !(w %in% SEVERITY_WAVES)) next
    tag <- paste(oc, w, sep = "_")
    # RESUME=1 skips models already fitted (an .rds in output/models), so an
    # interrupted run does not redo hours of work.
    if (RESUME && file.exists(file.path(DIR$models, paste0(tag, ".rds")))) next
    cost <- if (oc == "incidence") INC_PLAN$n_iter[INC_PLAN$wave == w] else
      CFG$n_iter_severity
    jobs[[length(jobs) + 1]] <- list(oc = oc, w = w, tag = tag,
                                     cost = as.numeric(cost))
  }
}
if (!length(jobs)) {
  log_msg("run_parallel: nothing to do, every model is already fitted")
  quit(save = "no")
}
jobs <- jobs[order(-vapply(jobs, `[[`, numeric(1), "cost"))]

loads <- rep(0, N); assign <- vector("list", N)
for (j in jobs) {
  i <- which.min(loads)
  assign[[i]] <- c(assign[[i]], list(j))
  loads[i] <- loads[i] + j$cost
}

log_msg("=== run_parallel: ", length(jobs), " models across ", N,
        " workers; estimated load (M iterations): ",
        paste(round(loads / 1e6, 1), collapse = " / "), " ===")

outdir <- file.path(ROOT, "output", "tables", "_workers")
if (RESUME && dir.exists(outdir)) {
  # keep the partial files of models already fitted; only clear the run marker
  # and the previous queues
  unlink(file.path(outdir, "RUNNING"))
  unlink(list.files(outdir, pattern = "^jobs\\.txt$", recursive = TRUE,
                    full.names = TRUE))
} else {
  unlink(outdir, recursive = TRUE)
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

specs <- vapply(seq_len(N), function(i) {
  paste(vapply(assign[[i]], function(j) paste(j$oc, j$w, sep = ":"),
               character(1)), collapse = ",")
}, character(1))

# One jobs.txt per worker so the progress panel can show the plan before the
# run starts and the current model while it runs.
for (i in seq_len(N)) {
  if (!nzchar(specs[i])) next
  wdir <- file.path(outdir, paste0("w", i))
  dir.create(wdir, recursive = TRUE, showWarnings = FALSE)
  writeLines(specs[i], file.path(wdir, "jobs.txt"))
}

# A failed model still leaves a .done marker, so on a resume run the marker must
# be cleared for every job about to be retried; otherwise the workers skip it and
# the wait loop returns at once.
for (tg in vapply(jobs, `[[`, character(1), "tag")) {
  old <- list.files(outdir, pattern = paste0("^", tg, "\\.done$"),
                    recursive = TRUE, full.names = TRUE)
  if (length(old)) unlink(old)
}

if (identical(Sys.getenv("PLAN_ONLY"), "1")) {
  for (i in seq_len(N)) if (nzchar(specs[i])) cat("w", i, ": ", specs[i], "\n", sep = "")
  log_msg("plan written; not launching (PLAN_ONLY=1)")
  quit(save = "no")
}

# ---- launch the workers ----------------------------------------------------
worker <- file.path(ROOT, "tools", "model_worker.R")
rscript <- file.path(R.home("bin"), "Rscript")
file.create(file.path(outdir, "RUNNING"))
procs <- list()

for (i in seq_len(N)) {
  spec <- specs[i]
  if (!nzchar(spec)) next
  wdir <- file.path(outdir, paste0("w", i))
  log_msg("worker ", i, ": ", spec)
  # parameters as arguments, not environment variables: system2(wait = FALSE)
  # does not pass `env` on Windows
  procs[[i]] <- system2(rscript,
                        c(shQuote(worker), shQuote(proj), shQuote(spec),
                          shQuote(wdir), shQuote(RUN_TAG_ARG),
                          shQuote(HACIN_SRC_ARG)),
                        wait = FALSE, stdout = FALSE, stderr = FALSE)
}

# ---- wait ------------------------------------------------------------------
expected <- vapply(jobs, `[[`, character(1), "tag")
repeat {
  done <- sub("\\.done$", "",
              list.files(outdir, pattern = "\\.done$", recursive = TRUE))
  if (all(expected %in% basename(done))) break
  Sys.sleep(20)
}
log_msg("all workers finished")
unlink(file.path(outdir, "RUNNING"))

# ---- merge -----------------------------------------------------------------
read_all <- function(suffix) {
  files <- list.files(outdir, pattern = paste0(suffix, "$"), recursive = TRUE,
                      full.names = TRUE)
  if (!length(files)) return(NULL)
  bind_rows(lapply(files, function(f) {
    utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  }))
}

# Rebuild the three tables from the saved model objects rather than from the
# partial files: after a resume the partials may hold an older fit of the same
# model (different covariates, different version), and merging them duplicates
# rows. The .rds in output/models is exactly one object per fitted model.
# Only the primary models are merged. Sensitivity scripts also save into
# output/models (14 saves <outcome>_Delta_vaccination.rds), and reading every
# .rds there duplicated the Delta coefficients.
expected_tags <- character(0)
for (oc in OUTCOMES) for (w in WAVES) {
  if (oc != "incidence" && !(w %in% SEVERITY_WAVES)) next
  expected_tags <- c(expected_tags, paste(oc, w, sep = "_"))
}
rds_files <- file.path(DIR$models, paste0(expected_tags, ".rds"))
rds_files <- rds_files[file.exists(rds_files)]
objs <- lapply(rds_files, readRDS)

coef_all <- bind_rows(lapply(objs, `[[`, "coefficients"))
area_all <- bind_rows(lapply(objs, function(o)
  cbind(tibble(outcome = o$outcome, wave = WAVE_LABELS[[o$wave]]), o$areas)))

# models that were not estimated keep their diagnostic row from the fit loop
diag_all <- bind_rows(lapply(objs, `[[`, "diagnostics"))
if (!is.null(diag_all) && nrow(diag_all)) {
  diag_all <- diag_all %>% distinct(outcome, wave, .keep_all = TRUE)
}

if (!is.null(coef_all)) write_csv_audit(coef_all, "coefficients.csv")
if (!is.null(diag_all)) write_csv_audit(diag_all, "diagnostics.csv")
if (!is.null(area_all)) write_csv_audit(area_all, "area_rr_summary.csv")

if (!is.null(diag_all)) {
  n_ok  <- sum(diag_all$status == "ok", na.rm = TRUE)
  n_net <- sum(grepl("^not estimated", diag_all$status))
  n_fail <- nrow(diag_all) - n_ok - n_net
  log_msg("run_parallel complete: ", n_ok, " fitted, ", n_net,
          " not estimated, ", n_fail, " failed")
  print(diag_all %>% mutate(across(c(rhat_max, ess_min), ~round(.x, 4))) %>%
          select(outcome, wave, rhat_max, ess_min, status))
}
