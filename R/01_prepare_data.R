# ===========================================================================
# 01_prepare_data.R -- build the single analysis dataset from the raw inputs.
#
# Inputs : data/raw/habana_panel_corregido.xlsx  (corrected panel; see docs/CORRECTIONS.md)
#          data/raw/shapefiles/AS_vuln_vs_casos_8_periodos.shp
#          data/raw/shapefiles/otras variables/AS_con_enlace_lista.shp
# Outputs: data/derived/analysis_dataset.rds  (sf, 574 rows)
#          data/derived/analysis_dataset.csv  (no geometry, for audit)
#          output/tables/data_checks.csv
#
# Reproducibility notes
#   - The legacy script joined the auxiliary shapefile by ROW POSITION. The
#     two files agree on Etiqueta (82/82) but disagree on one area NAME
#     ("Rafael Zilueta" vs "ROBERTO M ZULUETA"). The join here is by Etiqueta.
#   - The shapefile has no .prj (CRS undefined). Coordinates are metric
#     (UTM-like): total area = ~718 km2, consistent with Havana. Area is
#     therefore computed in native units / 1e6. See docs/DATA_DICTIONARY.md.
# ===========================================================================


source("config.R")
source(file.path(ROOT, "R", "00_setup.R"))

log_msg("--- 01_prepare_data ---")
invisible(lapply(FILES, assert_file))

# ---- 1. panel ---------------------------------------------------------------
panel_raw <- suppressMessages(
  read_excel(FILES$panel, sheet = 1, .name_repair = "minimal")
)

# The corrected workbook carries Prop_cancer twice: the original column and the
# one recomputed as counts / population * 100. `dat[["Prop_cancer"]]` would pick
# the first, i.e. the superseded column, so keep the last occurrence explicitly.
dup_cols <- unique(names(panel_raw)[duplicated(names(panel_raw))])
if (length(dup_cols)) {
  log_msg("duplicated columns, keeping last occurrence: ",
          paste(dup_cols, collapse = ", "))
  panel_raw <- panel_raw[, !duplicated(names(panel_raw), fromLast = TRUE)]
}

panel <- panel_raw %>%
  filter(!is.na(covid_variant), !is.na(ASalud)) %>%
  # the panel carries an empty geometry column and duplicate name columns;
  # drop them so the sf geometry and ASalud survive the join intact
  select(-any_of(c("geometry", "Municipio"))) %>%
  select(-ASalud)

log_msg("panel: ", nrow(panel_raw), " raw rows -> ", nrow(panel), " used")

assert(all(WAVES %in% unique(panel$covid_variant)),
       "not all seven waves are present in the panel")
assert(nrow(panel) == 82 * length(WAVES), "panel is not 82 areas x 7 waves")
assert(!anyDuplicated(panel[, c("codigo_as", "covid_variant")]),
       "duplicated area x wave keys in the panel")

# ---- 2. shapefiles ---------------------------------------------------------
shp <- read_sf(FILES$shp_main)
aux <- read_sf(FILES$shp_aux)

assert(nrow(shp) == 82, "main shapefile is not 82 areas")
assert(nrow(aux) == 82, "auxiliary shapefile is not 82 areas")
assert(!anyDuplicated(shp$Etiqueta), "duplicated Etiqueta in main shapefile")
assert(setequal(shp$Etiqueta, aux$Etiqueta),
       "Etiqueta sets differ between shapefiles")

# join auxiliary variables by CODE, never by position
aux_small <- aux %>%
  st_drop_geometry() %>%
  select(
    Etiqueta,
    PcvivT  = Ave_PcVivT, Equip   = Ave_EquipD, GradE   = Ave_GradxE,
    PcAAgua = Ave_PcAAgu, PcvivB  = Ave_PcVivB, PcvivD  = Ave_PcVivD,
    Pcviv1  = Ave_PcVi_1
  )

shp2 <- shp %>% left_join(aux_small, by = "Etiqueta")
assert(nrow(shp2) == 82, "join with auxiliary shapefile changed the row count")

# ---- 3. long dataset -------------------------------------------------------
data_shp <- shp2 %>%
  select(ASalud, Etiqueta, geometry, Equip, GradE, PcAAgua, PcvivD) %>%
  left_join(panel, by = c("Etiqueta" = "codigo_as"))

assert(nrow(data_shp) == 574, "joined dataset is not 574 rows")
assert(sum(is.na(data_shp$covid_variant)) == 0, "unmatched areas after join")

# ---- 3b. population and covariate source --------------------------------
# DECIDED (author): the correct population is `poblacion`. The panel columns
# poblacion_V, DensPob and DensPob_C are attached to the wrong areas and were
# removed from the corrected panel; poblacion_datoVac is kept for the audit
# only. `poblacion` is the offset for every model, the denominator of all
# rates, and (with the area computed from the geometry) the density.
assert(sum(is.na(data_shp$poblacion)) == 0, "population has missing values")
pop_total <- sum(data_shp$poblacion) / length(WAVES)
assert(pop_total > 1.9e6 && pop_total < 2.5e6,
       sprintf("implausible population total (%.0f)", pop_total))
log_msg("population total (poblacion): ", round(pop_total))

data_shp <- data_shp %>% mutate(poblacion_original = poblacion)

# Crowding: single definition (panel/census), decided by the author.
# The alternative derivation from the reviewed file was removed to avoid confusion.
assert("hacin" %in% names(data_shp), "crowding column missing from the panel")
log_msg("crowding: panel definition")

# ---- 4. derived variables --------------------------------------------------
# area in km2 (coordinates are metres); density as in the legacy analysis
data_shp <- data_shp %>%
  mutate(
    area_km2      = as.numeric(st_area(geometry)) / 1e6,
    density       = poblacion / area_km2,
    incidence_rate = cases_per_variant / poblacion * 1e5,   # per 100,000 (manuscript)
    mortality_rate = deaths_per_variant / poblacion * 1e5,
    incidence_cov  = log(cases_per_variant / poblacion * 1e5 + 1),   # covariate of the mortality model
    lethality      = deaths_per_variant / cases_per_variant * 100,
    cobertura_vac  = vacunados / poblacion * 100
  )

# ---- 4b. replication of the submitted analysis ----------------------------
# Uses the columns identified from Figure S1 of the submission. Vaccination
# coverage was computed on `poblacion` (before the offset is replaced).
if (SPEC == "submitted") {
  assert(all(c("DensPob_C", "CantDBprec_N", "CantDBprec", OFFSET_VAR) %in% names(data_shp)),
         "original panel columns needed for the replication are missing")
  data_shp <- data_shp %>%
    mutate(
      density    = DensPob_C,
      CantDBprec = if (SHANTY_VAR == "count") CantDBprec_N else CantDBprec,
      poblacion  = .data[[OFFSET_VAR]]
    )
  log_msg("SPEC = submitted; offset = ", OFFSET_VAR, "; shanty towns = ", SHANTY_VAR)
}

total_area <- sum(as.numeric(st_area(data_shp)) / 1e6) / length(WAVES)
assert(total_area > 600 && total_area < 900,
       sprintf("implausible total area (%.1f km2): check projection", total_area))

# ---- 5. checks for the audit trail ----------------------------------------
checks <- list(
  tibble(check = "rows", value = nrow(data_shp)),
  tibble(check = "areas", value = n_distinct(data_shp$ASalud)),
  tibble(check = "waves", value = n_distinct(data_shp$covid_variant)),
  tibble(check = "total_area_km2", value = round(total_area, 1)),
  tibble(check = "min_cases", value = min(data_shp$cases_per_variant)),
  tibble(check = "zero_case_rows", value = sum(data_shp$cases_per_variant == 0)),
  tibble(check = "zero_death_rows", value = sum(data_shp$deaths_per_variant == 0)),
  tibble(check = "vaccination_non_NA", value = sum(!is.na(data_shp$cobertura_vac))),
  tibble(check = "duplicated_keys",
         value = sum(duplicated(data_shp[, c("Etiqueta", "covid_variant")])))
) %>% bind_rows()

# zero-event counts by wave and outcome (Reviewer 1, point 1.15)
zero_events <- data_shp %>%
  st_drop_geometry() %>%
  group_by(wave = covid_variant) %>%
  summarise(
    areas            = n(),
    zero_cases       = sum(cases_per_variant == 0),
    zero_deaths      = sum(deaths_per_variant == 0),
    deaths_total     = sum(deaths_per_variant),
    cases_total      = sum(cases_per_variant),
    .groups = "drop"
  ) %>%
  mutate(wave = factor(wave, levels = WAVES, labels = WAVE_LABELS[WAVES])) %>%
  arrange(wave)

# ---- 6. save ---------------------------------------------------------------
save_rds(data_shp, "analysis_dataset.rds")
write_csv_audit(st_drop_geometry(data_shp), "analysis_dataset.csv", DIR$derived)
write_csv_audit(checks, "data_checks.csv")
write_csv_audit(zero_events, "zero_event_counts_by_wave.csv")

# input checksums
manifest <- tibble(
  file = unlist(FILES, use.names = TRUE),
  md5  = vapply(unlist(FILES), file_md5, character(1), USE.NAMES = FALSE)
)
write_csv_audit(manifest, "input_checksums.csv", DIR$derived)

log_msg("01_prepare_data complete: ", nrow(data_shp), " rows, ",
        n_distinct(data_shp$ASalud), " areas, ", n_distinct(data_shp$covid_variant), " waves")
print(zero_events)
