# ===========================================================================
# config.R -- single source of truth for paths, constants and options.
# Sourced by every script. Do not hard-code paths anywhere else.
# ===========================================================================

.this_file <- function() {
  # find the source() frame whose file is config.R itself, so that config.R
  # resolves the project root correctly even when sourced from a script in R/
  idx <- seq_len(sys.nframe())
  for (i in rev(idx)) {
    of <- sys.frame(i)$ofile
    if (!is.null(of) && basename(of) == "config.R") {
      return(normalizePath(of, mustWork = FALSE))
    }
  }
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[1]), mustWork = FALSE))
  normalizePath(".", mustWork = FALSE)
}

ROOT <- dirname(.this_file())

# RUN_TAG=<name> sends every derived file and output to a separate sub-folder,
# so sensitivity runs (e.g. the alternative crowding definition) never overwrite
# the main results. Empty = main run.
# SPEC selects the specification:
#   "corrected" (default): corrected panel and the decisions in docs/CORRECTIONS.md
#   "submitted"          : REPLICATION of the submitted analysis. Reads the original
#                          panel with the columns the manuscript used (identified from
#                          Figure S1: density = DensPob_C, crowding = hacin, ageing =
#                          ageing_index_cub of the panel; vaccination = vacunados /
#                          poblacion), the 9 covariates of Table 1, incidence only.
#                          The offset used originally is unknown, hence OFFSET_VAR:
#                          poblacion_datoVac | poblacion | poblacion_V
#                          and SHANTY_VAR: count (CantDBprec_N) | density (CantDBprec).
SPEC       <- Sys.getenv("SPEC", "corrected")
OFFSET_VAR <- Sys.getenv("OFFSET_VAR", "poblacion_datoVac")
SHANTY_VAR <- Sys.getenv("SHANTY_VAR", "count")
stopifnot(SPEC %in% c("corrected", "submitted"),
          OFFSET_VAR %in% c("poblacion_datoVac", "poblacion", "poblacion_V"),
          SHANTY_VAR %in% c("count", "density"))

# RUN_TAG=<name> (see below); a replication run is tagged automatically
RUN_TAG <- Sys.getenv("RUN_TAG")
if (SPEC == "submitted" && !nzchar(RUN_TAG)) {
  RUN_TAG <- paste0("submitted_", sub("poblacion_?", "pob_", OFFSET_VAR), "_shanty_", SHANTY_VAR)
}
.tagged <- function(...) {
  p <- file.path(ROOT, ...)
  if (nzchar(RUN_TAG)) file.path(p, RUN_TAG) else p
}

DIR <- list(
  raw     = file.path(ROOT, "data", "raw"),
  derived = .tagged("data", "derived"),
  tables  = .tagged("output", "tables"),
  figures = .tagged("output", "figures"),
  models  = .tagged("output", "models"),
  logs    = .tagged("output", "logs")
)
invisible(lapply(DIR, dir.create, recursive = TRUE, showWarnings = FALSE))

# --- input files ------------------------------------------------------------
# The restricted surveillance data are not distributed with this repository.
# If the corrected panel is absent from data/raw, the pipeline falls back to the
# synthetic dataset in data/synthetic (same files and structure, simulated
# values). Results produced in that mode are illustrative only.
PANEL_NAME <- "habana_panel_corregido.xlsx"
if (file.exists(file.path(DIR$raw, PANEL_NAME))) {
  DATA_DIR <- DIR$raw
} else {
  DATA_DIR <- file.path(ROOT, "data", "synthetic")
  message("Real data not found in data/raw: using the synthetic dataset in data/synthetic.")
  message("Results produced in this mode are illustrative and do not describe Havana.")
}

FILES <- list(
  # Corrected panel (see docs/CORRECTIONS.md). The original file
  # habana_covid_integration_by_variant_Wave_310825.xlsx is kept in data/raw
  # for the audit but is NOT used: several of its covariate columns are
  # attached to the wrong health area.
  panel    = if (SPEC == "submitted")
               file.path(DIR$raw, "habana_covid_integration_by_variant_Wave_310825.xlsx")
             else file.path(DATA_DIR, PANEL_NAME),
  shp_main = file.path(DATA_DIR, "shapefiles", "AS_vuln_vs_casos_8_periodos.shp"),
  shp_aux  = file.path(DATA_DIR, "shapefiles", "otras variables", "AS_con_enlace_lista.shp"),
  daily    = file.path(DATA_DIR, "covid_cases_havana_until_311222.xlsx"),
  serie    = file.path(DATA_DIR, "serie.xlsx")
)

# --- analysis options -------------------------------------------------------
CFG <- list(
  n_chains = 3L,
  # Chain length is set from effective sample size, not convention.
  # Incidence models carry 12 covariates plus a CAR effect on 82 areas and mix
  # slowly, so they get a long chain (thin 5 keeps the stored draws to 200,000
  # per chain). Severity models converge quickly and use a shorter chain.
  n_iter_incidence = 1000000L,
  burnin_incidence = 50000L,
  # thin 20 keeps ~47,500 stored draws per chain instead of 190,000: the
  # stored footprint of a model drops from ~2.3 GB to under 0.5 GB, which
  # matters on a 15.5 GB machine. Thinning does not change ESS.
  thin_incidence   = 20L,
  n_iter_severity  = 150000L,
  burnin_severity  = 15000L,
  thin_severity    = 1L,
  seed     = 12345L,

  # Covariate scaling. DECIDED: min-max to [0,1], as in the analysis code, so
  # every RR is per full observed range of the covariate. This keeps the
  # resubmission on the same scale as the submitted manuscript, so that any
  # difference between the two reflects the model, not the scaling.
  # NOTE: scaling changes the magnitude and the width of an RR interval but
  # never its statistical significance; see README.
  scaling  = "minmax",

  # Null-interval criteria. Primary is a 95% CrI excluding 1.
  null_primary = c(1.0, 1.0),
  null_band    = c(0.9, 1.1),
  null_sym     = c(0.8, 1.25),

  # Hotspot definition
  hotspot_rr        = 1.5,
  hotspot_min_waves = 5L,
  hotspot_prob      = 0.80
)

# QUICK_RUN=1 shortens the chains so the pipeline can be smoke-tested.
if (identical(Sys.getenv("QUICK_RUN"), "1")) {
  CFG$n_iter_incidence <- 20000L; CFG$burnin_incidence <- 5000L; CFG$thin_incidence <- 2L
  CFG$n_iter_severity  <- 20000L; CFG$burnin_severity  <- 5000L; CFG$thin_severity  <- 2L
  message("QUICK_RUN=1: short MCMC chains, results are NOT for publication.")
}

# ITER_OVERRIDE=<n> sets the number of iterations directly (testing only).
if (nzchar(Sys.getenv("ITER_OVERRIDE"))) {
  n <- as.integer(Sys.getenv("ITER_OVERRIDE")); b <- max(1000L, n %/% 10L)
  CFG$n_iter_incidence <- n; CFG$burnin_incidence <- b; CFG$thin_incidence <- 1L
  CFG$n_iter_severity  <- n; CFG$burnin_severity  <- b; CFG$thin_severity  <- 1L
  message("ITER_OVERRIDE=", n, ": results are NOT for publication.")
}

# --- waves ------------------------------------------------------------------
WAVES <- c("Wuhan1", "Wuhan2", "D614G", "Beta", "Delta", "Omicron1", "Omicron2")
WAVE_LABELS <- c(
  Wuhan1 = "Early", Wuhan2 = "Wuhan", D614G = "D614G", Beta = "Beta",
  Delta = "Delta", Omicron1 = "Omicron BA.1", Omicron2 = "Omicron BA.2"
)
OUTCOMES <- c("incidence", "mortality", "cfr")

# --- covariate sets (these reproduce Table 1 and Tables S1/S2) --------------
# One structural set for all three outcomes. Twelve covariates, all with
# variance inflation factor below 5. Of the pair of collinear housing variables
# (Equip and PcvivD, r = 0.90, VIF 26 and 38) PcvivD is retained, because it is
# the housing-typology measure named in the Methods. CMF density enters as
# Dens_x_CMF, the version aligned with the area codes; Dens_x_CMF_A is not used.
# Pob60aPerS was dropped because it duplicates the ageing index.
COV_INCIDENCE <- c(
  "density", "GradE", "ageing_index_cub", "hacin", "CantDBprec",
  "Dis_x_CMF_A", "PobAt_CMF_A", "Aglomeraci", "casasAlq",
  "PcAAgua", "PcvivD", "Dens_x_CMF"
)
COV_VACCINE <- "cobertura_vac"   # Delta incidence model only

COV_COMORBID <- c(
  "Prop_Hip", "Prop_Diab", "Prop_Card", "Prop_Epoc", "Prop_cancer", "Prop_Asma"
)
# Severity drops population density. The twelve structural covariates plus the
# six comorbidities reach a variance inflation factor of 5.6-5.7, always on
# density: it is already 78% predictable from the other structural covariates
# and the six prevalence variables add the few points of R2 that push it over.
# Without density the maximum VIF falls to 3.6 in case fatality and 3.9 in
# mortality, so every model stays below the threshold stated in the Methods.
# Density remains in the incidence models, where the maximum is 4.5.
COV_SEVERITY <- c(setdiff(COV_INCIDENCE, "density"), COV_COMORBID)

# Reduced structural set, matching the covariate list printed in Table 1 of the
# submitted manuscript. Used as a sensitivity analysis.
COV_INCIDENCE_REDUCED <- c(
  "density", "GradE", "ageing_index_cub", "hacin", "CantDBprec",
  "Dis_x_CMF_A", "PobAt_CMF_A", "Aglomeraci", "casasAlq"
)

# Replication of the submitted analysis: 9 covariates of Table 1, incidence only
# (the covariates of the submitted severity models are not documented).
if (SPEC == "submitted") {
  COV_INCIDENCE <- COV_INCIDENCE_REDUCED
  COV_SEVERITY  <- COV_INCIDENCE_REDUCED
  OUTCOMES      <- "incidence"
}

# Severity (mortality and case fatality) models are fitted for the first five
# waves (Early, Wuhan, D614G, Beta, Delta). DECIDED (author): only the two
# Omicron waves are left out. The death count is no longer a filter: waves
# below MIN_DEATHS_SEVERITY deaths are still fitted but flagged "sparse
# outcome" in diagnostics.csv, and their estimates are imprecise (Early has 48
# deaths and Wuhan 20 for 18 covariates) and should be read as exploratory.
SEVERITY_WAVES <- c("Wuhan1", "Wuhan2", "D614G", "Beta", "Delta")
MIN_DEATHS_SEVERITY <- 100L

# Iterations per incidence wave, set from the effective sample sizes obtained in
# earlier runs: the early and Omicron waves mix well at 1M, D614G and Beta need
# 3M, and Delta needs 6M. Thinning is chosen so that the stored draws stay within
# memory; it does not affect ESS.
INC_PLAN <- data.frame(
  wave   = c("Wuhan1", "Wuhan2", "D614G", "Beta", "Delta", "Omicron1", "Omicron2"),
  n_iter = c(1000000, 1000000, 3000000, 3000000, 6000000, 1000000, 1000000),
  burnin = c(100000,  100000,  300000,  300000,  600000,  100000,  100000),
  thin   = c(10,       10,       20,       20,       40,       10,       10)
)

# Variables described in the Methods but excluded before fitting
# (collinearity / instability). Retained here only for the audit trail.
# Candidate variables not used in the final models.
#   Equip        collinear with PcvivD (r = 0.90, VIF 26 and 38); PcvivD retained
#                because it is the housing-typology measure named in the Methods
#   Dens_x_CMF_A the version in the panel is not aligned with the area codes
#                (r = 0.001 with distance to the nearest CMF, against -0.36 for
#                the aligned version Dens_x_CMF, which is used instead)
COV_EXCLUDED <- c("Equip", "Dens_x_CMF_A")

# --- display labels ---------------------------------------------------------
LABELS <- c(
  density = "Population density", GradE = "Education level",
  ageing_index_cub = "Ageing index", hacin = "Crowding",
  CantDBprec = "Shanty towns", Dis_x_CMF_A = "Distance to nearest CMF",
  PobAt_CMF_A = "Population served per CMF", Aglomeraci = "Agglomeration sites",
  casasAlq = "Rental housing density", cobertura_vac = "Vaccination coverage",
  PcAAgua = "Piped-water access", PcvivD = "Precarious housing",
  incidence_cov = "Cumulative incidence (log, same wave)",
  Dens_x_CMF = "CMF density", Prop_Hip = "Hypertension", Prop_Diab = "Diabetes", Prop_Card = "Cardiovascular disease",
  Prop_Epoc = "COPD", Prop_cancer = "Cancer", Prop_Asma = "Asthma"
)

# Wave-specific covariate: crude cumulative incidence of the SAME wave, as
# log(incidence per 100,000 + 1), standardised within each wave (03_standardise.R).
#   mortality: adjusted for incidence (default TRUE). Mortality = incidence x
#     case fatality, so without it the mortality coefficients mix transmission
#     and severity. With it they are effects on deaths given the observed
#     incidence (not mediated through incidence).
#   case fatality: NOT adjusted by default. Cases are already the offset, so
#     incidence (cases/population) would share the denominator (mathematical
#     coupling) and case under-ascertainment would create spurious
#     associations. Available as a sensitivity (20_severity_incidence_adjustment.R).
COV_WAVE_VARYING <- "incidence_cov"
CFG$mortality_adjust_incidence <- TRUE
CFG$cfr_adjust_incidence       <- FALSE

covariates_for <- function(outcome, wave = NULL) {
  if (outcome == "incidence") {
    v <- COV_INCIDENCE
    if (identical(wave, "Delta")) v <- c(v, COV_VACCINE)
    return(v)
  }
  v <- COV_SEVERITY
  if ((outcome == "mortality" && isTRUE(CFG$mortality_adjust_incidence)) ||
      (outcome == "cfr"       && isTRUE(CFG$cfr_adjust_incidence))) {
    v <- c(v, COV_WAVE_VARYING)
  }
  v
}

outcome_cols <- function(outcome) {
  switch(outcome,
    incidence = list(y = "cases_per_variant", offset = "poblacion"),
    mortality = list(y = "deaths_per_variant", offset = "poblacion"),
    cfr       = list(y = "deaths_per_variant", offset = "cases_per_variant")
  )
}

label_of <- function(x) ifelse(x %in% names(LABELS), LABELS[x], x)
