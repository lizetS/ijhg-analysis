#!/usr/bin/env Rscript
# ===========================================================================
# run_all.R -- rebuild every reported result from the raw inputs.
#
# Usage
#   Rscript run_all.R                      # main analysis only
#   Rscript run_all.R --with-sensitivity   # plus the reviewer-requested analyses
#
#   QUICK_RUN=1 Rscript run_all.R          # short chains, smoke test only
#
# The pipeline is deterministic: one seed, fixed inputs, no manual steps.
# ===========================================================================

# ---- make the working directory the repository folder (this file's folder) ----
# Works with Rscript, with "Source" in RStudio and when sourced from the console.
.script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) return(dirname(normalizePath(f[1])))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(dirname(normalizePath(of)))
  }
  getwd()
}
setwd(.script_dir())
if (!file.exists("config.R"))
  stop("config.R not found in ", getwd(), ". Set the working directory to the ",
       "repository folder (the one containing config.R and run_all.R), e.g. ",
       "setwd('C:/path/to/reproducible_analysis')")

args <- commandArgs(trailingOnly = TRUE)
with_sens <- "--with-sensitivity" %in% args

t_start <- Sys.time()
source("config.R")
source(file.path(ROOT, "R", "00_setup.R"))

# record the environment in which the results were produced
writeLines(capture.output(sessionInfo()),
           file.path(DIR$logs, "session_info.txt"))

scripts <- c(
  "R/01_prepare_data.R",
  "R/02_waves_descriptives.R",
  "R/03_standardise.R",
  "R/04_fit_models.R",
  "R/05_tables.R",
  "R/21_collinearity.R",
  "R/23_table2_covariates.R",
  "R/06_null_interval.R",
  "R/07_hotspots.R",
  "R/09_figures.R",
  "R/16_figures_descriptive.R",
  "R/19_figure_S1.R",
  "R/17_figure1.R",
  "R/25_extra_outputs.R"
)
if (with_sens) scripts <- append(scripts, c("R/08_sensitivity.R", "R/14_vaccination_severity.R", "R/18_regularised_severity.R", "R/20_severity_incidence_adjustment.R", "R/26_bym.R"), after = 7)

log_msg("=== run_all start ===")
log_msg("R ", getRversion(), "; scaling = ", CFG$scaling,
        "; chains = ", CFG$n_chains, "; iterations = ", CFG$n_iter)

for (s in scripts) {
  log_msg("> ", s)
  t0 <- Sys.time()
  source(file.path(ROOT, s))
  log_msg("< ", s, " (", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), " min)")
}

# manifest of everything produced
outs <- list.files(c(DIR$tables, DIR$figures, DIR$models), recursive = FALSE,
                   full.names = TRUE)
manifest <- data.frame(
  file = sub(paste0(ROOT, "/"), "", outs),
  size_bytes = file.info(outs)$size,
  md5 = unname(tools::md5sum(outs))
)
write.csv(manifest, file.path(DIR$logs, "output_manifest.csv"), row.names = FALSE)

log_msg("=== run_all complete in ",
        round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1),
        " minutes; ", nrow(manifest), " outputs ===")
