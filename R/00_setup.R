# ===========================================================================
# 00_setup.R -- packages and small helpers shared by all scripts.
# Sourced after config.R; must not source itself.
# ===========================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readxl)
  library(writexl)
  library(sf)
  library(spdep)
  library(scales)
  library(CARBayes)
  library(coda)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

# MASS (pulled in by CARBayes) masks dplyr::select and stats::filter.
select <- dplyr::select
filter <- dplyr::filter
count  <- dplyr::count

LOG_FILE <- file.path(DIR$logs, "run.log")

log_msg <- function(...) {
  txt <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...))
  cat(txt, "\n", sep = "")
  cat(txt, "\n", sep = "", file = LOG_FILE, append = TRUE)
}

assert <- function(cond, msg) {
  if (!isTRUE(cond)) stop("ASSERTION FAILED: ", msg, call. = FALSE)
  invisible(TRUE)
}

assert_file <- function(path) {
  assert(file.exists(path), paste("missing input:", path))
}

file_md5 <- function(path) unname(tools::md5sum(path))

# Write a data frame to CSV and register it in the output manifest.
write_csv_audit <- function(df, filename, dir = DIR$tables) {
  path <- file.path(dir, filename)
  utils::write.csv(df, path, row.names = FALSE, na = "")
  log_msg("wrote ", path)
  invisible(path)
}

save_rds <- function(obj, filename, dir = DIR$derived) {
  path <- file.path(dir, filename)
  saveRDS(obj, path)
  log_msg("wrote ", path)
  invisible(path)
}

read_rds <- function(filename, dir = DIR$derived) {
  path <- file.path(dir, filename)
  assert_file(path)
  readRDS(path)
}

# Build the model formula for one outcome and covariate set.
build_formula <- function(covs, outcome) {
  oc <- outcome_cols(outcome)
  rhs <- paste(c(sprintf("offset(log(%s))", oc$offset), covs), collapse = " + ")
  stats::as.formula(paste(oc$y, "~", rhs))
}
