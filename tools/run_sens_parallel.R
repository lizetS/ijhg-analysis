# ===========================================================================
# run_sens_parallel.R -- run sections A, A2, B, C and D of 08_sensitivity.R
# across several processes and rebuild their output tables.
#
# 08 as written is a single bucle: 61 fits, 424 minutes on one core. The jobs
# here are independent of each other, so they are spread over N workers and the
# tables are merged at the end.
#
#   N_WORKERS=12 Rscript tools/run_sens_parallel.R
#
# Outputs: TableS4_incidence_comorbid.csv, sensitivity_reduced_covariates.csv,
#          sensitivity_binomial_cfr.csv, sensitivity_weights.csv,
#          sensitivity_priors.csv
# ===========================================================================

proj <- Sys.getenv("PROJ", unset = getwd())
setwd(proj)
source("config.R")
source(file.path(ROOT, "R", "00_setup.R"))

N <- as.integer(Sys.getenv("N_WORKERS", "8"))
stopifnot(!is.na(N), N >= 1)

log_msg("=== run_sens_parallel: ", N, " workers ===")

# ---- job list ---------------------------------------------------------------
# cost is in iterations, taken from the same plan the primary models use
cost_of <- function(oc, wv) if (oc == "incidence")
  as.numeric(INC_PLAN$n_iter[INC_PLAN$wave == wv]) else CFG$n_iter_severity

jobs <- list()
add <- function(sec, oc, waves, variant = "-", factor = 1) {
  for (wv in waves) {
    jobs[[length(jobs) + 1]] <<- list(
      spec = paste(sec, oc, wv, variant, sep = ";"),
      tag  = paste(sec, oc, wv, variant, sep = "_"),
      cost = cost_of(oc, wv) * factor)
  }
}

# SECTIONS=A,C lets a partial rerun touch only those tables
SECTIONS <- strsplit(Sys.getenv("SECTIONS", "A,A2,B,C,D"), ",", fixed = TRUE)[[1]]
want <- function(s) s %in% SECTIONS

if (want("A"))  add("A",  "incidence", WAVES)                    # comorbidities
if (want("A2")) add("A2", "incidence", WAVES, factor = 0.6)      # reduced set
if (want("B"))  add("B",  "cfr", SEVERITY_WAVES, variant = "-", factor = 0.15)
if (want("C")) {
  add("C", "incidence", WAVES, "queen");  add("C", "incidence", WAVES, "rook")
  add("C", "incidence", WAVES, "queen2"); add("C", "incidence", WAVES, "knn4")
}
if (want("D")) {
  add("D", "incidence", WAVES, "default")
  add("D", "incidence", WAVES, "ig_half")
  add("D", "incidence", WAVES, "ig_vague")
}

jobs <- jobs[order(-vapply(jobs, `[[`, numeric(1), "cost"))]
cat("trabajos:", length(jobs), "\n")

loads <- rep(0, N); assign <- vector("list", N)
for (j in jobs) {
  i <- which.min(loads)
  assign[[i]] <- c(assign[[i]], list(j))
  loads[i] <- loads[i] + j$cost
}
log_msg("carga estimada por worker (M iteraciones): ",
        paste(round(loads / 1e6, 1), collapse = " / "))

outdir <- file.path(ROOT, "output", "tables", "_sens")

if (identical(Sys.getenv("MERGE_ONLY"), "1")) {
  log_msg("MERGE_ONLY=1: reusing the outputs already in ", outdir)
} else {
unlink(outdir, recursive = TRUE); dir.create(outdir, recursive = TRUE)
# ---- launch -----------------------------------------------------------------
worker <- file.path(ROOT, "tools", "sens_worker.R")
rscript <- file.path(R.home("bin"), "Rscript")
for (i in seq_len(N)) {
  if (!length(assign[[i]])) next
  spec <- paste(vapply(assign[[i]], `[[`, character(1), "spec"), collapse = ",")
  wdir <- file.path(outdir, paste0("w", i)); dir.create(wdir, recursive = TRUE)
  system2(rscript, c(shQuote(worker), shQuote(proj), shQuote(spec), shQuote(wdir)),
          wait = FALSE, stdout = FALSE, stderr = FALSE)
  log_msg("worker ", i, ": ", length(assign[[i]]), " trabajos")
}

expected <- vapply(jobs, `[[`, character(1), "tag")
repeat {
  done <- sub("\\.done$", "", list.files(outdir, pattern = "\\.done$", recursive = TRUE))
  if (all(expected %in% basename(done))) break
  Sys.sleep(30)
}
log_msg("todos los workers terminaron")
}

# ---- merge ------------------------------------------------------------------
prim <- utils::read.csv(file.path(DIR$tables, "coefficients.csv"),
                        stringsAsFactors = FALSE)

read_tag <- function(tag) {
  f <- list.files(outdir, pattern = paste0("^", tag, "\\.csv$"), recursive = TRUE,
                  full.names = TRUE)
  if (!length(f)) return(NULL)
  utils::read.csv(f[1], stringsAsFactors = FALSE, check.names = FALSE)
}
all_jobs <- lapply(jobs, function(j) {
  d <- read_tag(j$tag); if (is.null(d)) return(NULL)
  f <- strsplit(j$spec, ";", fixed = TRUE)[[1]]
  data.frame(d, section = f[1], variant = f[4], stringsAsFactors = FALSE)
})

grab <- function(sec, var = NULL) {
  idx <- vapply(jobs, function(j) {
    f <- strsplit(j$spec, ";", fixed = TRUE)[[1]]
    f[1] == sec && (is.null(var) || f[4] == var)
  }, logical(1))
  bind_rows(all_jobs[idx])
}

# A: incidence with the six comorbidities
ga <- grab("A")
if (!is.null(ga) && nrow(ga)) write_csv_audit(ga %>% select(outcome, wave, variable, rr, lower, upper, pr_gt_1),
                              "TableS4_incidence_comorbid.csv")

# B: binomial case fatality, side by side with the Poisson primary model
prim_cfr <- prim %>% filter(outcome == "cfr") %>% mutate(likelihood = "poisson")
bin <- grab("B")
if (!is.null(bin) && nrow(bin)) {
  cmp_b <- bind_rows(prim_cfr, bin %>% mutate(likelihood = "binomial")) %>%
    filter(variable != "(Intercept)") %>%
    transmute(wave, variable, likelihood,
              rr = round(rr, 3), lower = round(lower, 3), upper = round(upper, 3),
              excludes_1 = lower > 1 | upper < 1) %>%
    tidyr::pivot_wider(names_from = likelihood,
                       values_from = c(rr, lower, upper, excludes_1))
  write_csv_audit(cmp_b, "sensitivity_binomial_cfr.csv")
}

# C and D keep their variant column, as 08 wrote them
gc_ <- grab("C")
if (!is.null(gc_) && nrow(gc_)) write_csv_audit(gc_ %>% rename(weights = variant) %>%
                  select(weights, outcome, wave, variable, rr, lower, upper, pr_gt_1),
                "sensitivity_weights.csv")
gd <- grab("D")
if (!is.null(gd) && nrow(gd)) write_csv_audit(gd %>% rename(prior = variant) %>%
                  select(prior, outcome, wave, variable, rr, lower, upper, pr_gt_1),
                "sensitivity_priors.csv")

# A2: comparison against the primary models
red <- grab("A2")
if (nrow(red)) {
  cmp <- bind_rows(prim %>% filter(outcome == "incidence") %>% mutate(model = "full"),
                   red %>% mutate(model = "reduced")) %>%
    select(model, wave, variable, rr, lower, upper) %>%
    tidyr::pivot_wider(names_from = model, values_from = c(rr, lower, upper))
  write_csv_audit(cmp, "sensitivity_reduced_covariates.csv")
}

log_msg("run_sens_parallel complete")
