# Data dictionary

Two raw sources are combined:

1. **Panel** — `data/raw/habana_covid_integration_by_variant_Wave_310825.xlsx`
   (one row per health area × wave; 574 usable rows after removing 1,048,001 empty rows).
2. **Geometry and housing variables** — `shapefiles/AS_vuln_vs_casos_8_periodos.shp`
   and `shapefiles/otras variables/AS_con_enlace_lista.shp`, joined on `Etiqueta`.

## Outcomes and offsets

| Variable | Definition | Source |
|---|---|---|
| `cases_per_variant` | confirmed cases, area × wave | MINSAP panel |
| `deaths_per_variant` | deaths, area × wave | MINSAP panel |
| `poblacion` | **corrected population**, area × wave | `poblacion_datoVac` in the panel |
| `poblacion_original` | uncorrected population, kept for audit | `poblacion` in the panel |

Model offsets: `log(poblacion)` for incidence and mortality; `log(cases_per_variant)` for
case fatality. No area has zero cases, so no row is lost to a zero offset.

## Derived variables

| Variable | Definition |
|---|---|
| `area_km2` | polygon area / 1e6 (coordinates are metres; CRS undefined) |
| `density` | `poblacion / area_km2`, using the corrected population |
| `incidence_rate` | `cases / poblacion * 1e5` (descriptive only) |
| `mortality_rate` | `deaths / poblacion * 1e5` (descriptive only) |
| `lethality` | `deaths / cases * 100` (descriptive only) |
| `cobertura_vac` | `vacunados / poblacion * 100`; defined in the Delta wave only |

## Note on the distance variable

Two distance variables exist in the panel: `Dis_x_CMF_A` (used here, and the one behind the
published Table 1) and `Dist_x_AS` (used by the July 2025 scripts). They are different
measurements: raw correlation 0.64, with large per-area discrepancies (e.g. Managua 481 vs
1,585). **`Dis_x_CMF_A` is the correct measure** (confirmed by the author). Its values are in
metres; because covariates are min–max rescaled, the unit does not affect any reported RR,
but the Methods must not describe it in kilometres.

## Covariates

All are standardised across the 82 areas (`R/03_standardise.R`), so every reported RR is
per one-standard-deviation increment when `CFG$scaling = "z"`.

| Code | Label | Definition | Source |
|---|---|---|---|
| `density` | Population density | inhabitants / km² | 2012 census |
| `GradE` | Education level | mean % secondary or higher (`Ave_GradxE`) | 2012 census |
| `ageing_index_cub` | Ageing index | % aged ≥ 60 | 2012 census |
| `hacin` | Crowding | residents per bedroom | 2012 census |
| `CantDBprec` | Shanty towns | count of informal settlements | 2012 census |
| `Dis_x_CMF_A` | Distance to nearest CMF | metres (100–1,726) | MINSAP 2020 |
| `PobAt_CMF_A` | Population served per CMF | residents per physician | MINSAP 2020 |
| `Aglomeraci` | Agglomeration sites | count of commercial/service facilities | 2012 census |
| `Dens_x_CMF` | CMF density | family medical centres per km² (aligned version) | MINSAP 2020 |
| `casasAlq` | Rental housing density | hotels + tourist rentals per km² | ONAT 2020 |
| `cobertura_vac` | Vaccination coverage | % completed primary schedule | ANDARIEGO (Delta only) |
| `Prop_Hip` | Hypertension | registered per 1,000 | primary care registry 2020 |
| `Prop_Diab` | Diabetes | registered per 1,000 | primary care registry 2020 |
| `Prop_Card` | Cardiovascular disease | registered per 1,000 | primary care registry 2020 |
| `Prop_Epoc` | COPD | registered per 1,000 | primary care registry 2020 |
| `Prop_cancer` | Cancer | registered per 1,000 | primary care registry 2020 |
| `Prop_Asma` | Asthma | registered per 1,000 | primary care registry 2020 |

## Variables described but not modelled

`Equip` (equipment-to-resident ratio) and `Dens_x_CMF_A` (CMF density, unaligned version).

`Equip` is collinear with `PcvivD` (r = 0.90; VIF 26 and 38 when both are entered); `PcvivD`
is retained because it is the housing-typology measure named in the Methods.
`Dens_x_CMF_A` is not aligned with the area codes — it is uncorrelated with distance to the
nearest CMF (r = 0.001), whereas the aligned `Dens_x_CMF` correlates at −0.36 — so the
aligned variable is used instead.

`PcAAgua` (piped water) and `PcvivD` (precarious housing) **are** in the final models.

## Wave coding

| Code | Label | Lineage |
|---|---|---|
| `Wuhan1` | Early | before variant characterisation |
| `Wuhan2` | Wuhan | ancestral lineage |
| `D614G` | D614G | D614G |
| `Beta` | Beta | B.1.351 |
| `Delta` | Delta | B.1.617.2 |
| `Omicron1` | Omicron BA.1 | BA.1 |
| `Omicron2` | Omicron BA.2 | BA.2 |

The code names are retained from the source panel for traceability; all outputs use the
manuscript labels.
