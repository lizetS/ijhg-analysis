#!/usr/bin/env Rscript
# ===========================================================================
# make_synthetic_data.R -- build the synthetic test dataset shipped with the
# repository.
#
# Nothing in it corresponds to a real health area or a real person: the
# geography is a Voronoi tessellation of a 30 x 24 km rectangle and every
# covariate and outcome is simulated. It has the same structure as the
# restricted inputs (82 areas x 7 waves, same files and columns) so that the
# whole pipeline can be run and tested.
#
# Usage: Rscript tools/make_synthetic_data.R
# ===========================================================================

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(tidyr); library(writexl)
})

a <- commandArgs(trailingOnly = FALSE)
f <- sub("^--file=", "", a[grepl("^--file=", a)])
ROOT <- if (length(f)) dirname(dirname(normalizePath(f[1]))) else getwd()
setwd(ROOT)
syn <- file.path("data", "synthetic")
dir.create(file.path(syn, "shapefiles", "otras variables"),
           recursive = TRUE, showWarnings = FALSE)
set.seed(20261101)

n <- 82L
WAVES <- c("Wuhan1", "Wuhan2", "D614G", "Beta", "Delta", "Omicron1", "Omicron2")
START <- as.Date(c("2020-03-16", "2020-07-05", "2020-11-06", "2021-03-27",
                   "2021-07-03", "2021-12-19", "2022-03-13"))
END <- c(START[-1] - 1, as.Date("2022-05-31"))
RATE <- c(60, 100, 1400, 2400, 3800, 240, 115)      # cases per 100k
CFR  <- c(3.7, 0.9, 0.57, 0.83, 0.60, 0.07, 0.08) / 100

# ---- geography: Voronoi tessellation of 30 x 24 km -------------------------
pts <- data.frame(x = runif(n, 0, 30000), y = runif(n, 0, 24000),
                  Etiqueta = seq_len(n), ASalud = sprintf("AREA %02d", seq_len(n)))
p <- st_as_sf(pts, coords = c("x", "y"))
env <- st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 30000, ymax = 24000)))
vor <- st_voronoi(st_union(p), envelope = env)
cells <- st_collection_extract(vor, "POLYGON")
polys <- st_intersection(st_sf(geometry = cells), env)  # clip to the 30 x 24 km rectangle
hit <- st_intersects(polys, p)
who <- vapply(hit, function(i) if (length(i)) i[1] else NA_integer_, integer(1))
polys$Etiqueta <- pts$Etiqueta[who]
polys$ASalud <- pts$ASalud[who]
polys <- polys[order(polys$Etiqueta), ]
stopifnot(nrow(polys) == n, !anyNA(polys$Etiqueta))
polys$Municipio <- rep(sprintf("MUNICIPIO %d", 1:15), length.out = n)
st_crs(polys) <- 32617   # UTM 17N, metres (the real file has no .prj)

dir.create(file.path(syn, "shapefiles"), recursive = TRUE, showWarnings = FALSE)
st_write(polys[, c("Etiqueta", "ASalud", "Municipio", "geometry")],
         file.path(syn, "shapefiles", "AS_vuln_vs_casos_8_periodos.shp"),
         delete_dsn = TRUE, quiet = TRUE)

# ---- auxiliary shapefile: census-type aggregates ---------------------------
aux <- polys %>% st_drop_geometry() %>%
  mutate(
    Ave_PcVivT = runif(n, 80, 100),
    Ave_PcVivB = runif(n, 0, 20),
    Ave_PcVivD = pmin(99.7, pmax(45, rnorm(n, 95, 8))),
    Ave_PcVi_1 = runif(n, 0, 30),
    Ave_EquipD = rgamma(n, 2, 0.10),
    Ave_GradxE = pmin(70, pmax(20, rnorm(n, 48, 8))),
    Ave_PcAAgu = pmin(50, rgamma(n, 3, 0.2))
  )
dir.create(file.path(syn, "shapefiles", "otras variables"), recursive = TRUE, showWarnings = FALSE)
st_write(st_sf(aux, geometry = polys$geometry),
         file.path(syn, "shapefiles", "otras variables", "AS_con_enlace_lista.shp"),
         delete_dsn = TRUE, quiet = TRUE)

# ---- area covariates -------------------------------------------------------
area_km2 <- as.numeric(st_area(polys)) / 1e6
pop <- rlnorm(n, log(2317226 / n), 0.35)
pop <- round(pop / sum(pop) * 2317226)
base <- tibble(
  codigo_as = polys$Etiqueta,
  ASalud = polys$ASalud,
  Municipio = polys$Municipio,
  poblacion = pop,
  area_km2 = area_km2,
  hacin = pmin(3.2, pmax(1.3, rnorm(n, 1.84, 0.30))),
  ageing_index_cub = pmin(40, pmax(15, rnorm(n, 27, 5.5))),
  CantDBprec = rpois(n, 0.9),
  Dis_x_CMF_A = pmin(1700, pmax(100, rgamma(n, 2, 0.005))),
  PobAt_CMF_A = pmin(18000, pmax(200, rlnorm(n, log(1800), 0.8))),
  Dens_x_CMF = pmin(80, pmax(0.5, rgamma(n, 1.2, 0.10))),
  Aglomeraci = pmax(0, rgamma(n, 2, 0.30)),
  casasAlq = pmax(0, rgamma(n, 0.6, 0.015)),
  Prop_Hip = pmin(40, pmax(0.3, rnorm(n, 16.6, 6.4))),
  Prop_Diab = pmin(40, pmax(0.3, rgamma(n, 1.2, 6))),
  Prop_Card = pmin(20, pmax(0.1, rgamma(n, 2, 2))),
  Prop_Epoc = pmin(5.5, pmax(0.05, rgamma(n, 1, 1))),
  Prop_cancer = pmin(40, pmax(0.8, rgamma(n, 1, 6.4))),
  Prop_Asma = pmin(35, pmax(0.05, rgamma(n, 1, 4)))
)

z <- function(x) as.numeric(scale(x))
# modest, wave-specific effects so that the models recover plausible estimates
mult_wave <- function(w) {
  eta <- 0.25 * z(base$hacin) + 0.20 * z(base$casasAlq) +
    0.30 * z(base$Prop_Diab) * (w %in% c("Omicron1", "Omicron2")) -
    0.30 * z(base$area_km2) * (w %in% c("Omicron1", "Omicron2")) +
    rnorm(n, 0, 0.15)
  exp(eta - mean(eta))
}

rows <- vector("list", length(WAVES))
for (i in seq_along(WAVES)) {
  m <- mult_wave(WAVES[i])
  lam <- base$poblacion * (RATE[i] / 1e5) * m
  cases <- rpois(n, lam)
  cases <- pmax(cases, 1L)
  p_death <- pmin(0.6, CFR[i] * exp(0.25 * z(base$Prop_Diab) +
                                     0.20 * z(base$Prop_Epoc) - 1))
  deaths <- rbinom(n, cases, p_death)
  rows[[i]] <- base %>%
    mutate(
      covid_variant = WAVES[i],
      start_date = START[i],
      end_date = END[i],
      cases_per_variant = cases,
      deaths_per_variant = deaths,
      vacunados = if (WAVES[i] == "Delta") round(poblacion * runif(n, 0, 0.70)) else NA_real_
    )
}
panel <- bind_rows(rows) %>%
  mutate(CantDBprec_N = CantDBprec) %>%
  select(ASalud, Municipio, covid_variant, codigo_as, poblacion, vacunados,
         cases_per_variant, deaths_per_variant, start_date, end_date,
         hacin, ageing_index_cub, CantDBprec, CantDBprec_N, Dis_x_CMF_A, PobAt_CMF_A,
         Dens_x_CMF, Aglomeraci, casasAlq,
         Prop_Hip, Prop_Diab, Prop_Card, Prop_Epoc, Prop_cancer, Prop_Asma)

notas <- tibble(
  columna = c("todas"),
  accion = c("SYNTHETIC TEST DATA: simulated values, not real Havana data")
)
write_xlsx(list(panel_corregido = panel, cambios = notas, pendientes = notas),
           file.path(syn, "habana_panel_corregido.xlsx"))

# ---- daily series and line list -------------------------------------------
dates <- seq(as.Date("2020-03-11"), as.Date("2022-05-31"), by = "day")
inc <- numeric(length(dates))
for (i in seq_along(WAVES)) {
  idx <- which(dates >= START[i] & dates <= END[i])
  k <- length(idx)
  shape <- sin(seq(0, pi, length.out = k))^2
  cases_w <- sum(panel$cases_per_variant[panel$covid_variant == WAVES[i]])
  inc[idx] <- cases_w * shape / sum(shape)
}
inc <- rpois(length(dates), inc)
wave_of_day <- pmax(1L, findInterval(dates, START))
dth <- rpois(length(dates), inc * CFR[wave_of_day])
dth <- pmin(dth, inc)
serie <- tibble(Fecha = dates, Incidence = inc, tiempo = seq_along(dates),
                deaths = dth,
                Mortality = dth / 2317226 * 1e5,
                CFR = ifelse(inc > 0, dth / inc * 100, 0))
write_xlsx(serie, file.path(syn, "serie.xlsx"))

n_ll <- 40000L
pick <- sample(seq_along(dates), n_ll, replace = TRUE, prob = pmax(inc, 0.01))
linelist <- tibble(
  date_conf = dates[pick],
  Health_Area = sample(base$ASalud, n_ll, replace = TRUE),
  CODEAREA = sample(base$codigo_as, n_ll, replace = TRUE),
  evolucion = ifelse(runif(n_ll) < 0.01, "FALLECIDO", NA_character_),
  Date_death = as.Date(NA)
)
write_xlsx(linelist, file.path(syn, "covid_cases_havana_until_311222.xlsx"))

cat("synthetic dataset written to", syn, "\n")
cat("areas:", n, " waves:", length(WAVES), " panel rows:", nrow(panel), "\n")
cat("total cases:", sum(panel$cases_per_variant), " deaths:", sum(panel$deaths_per_variant), "\n")
