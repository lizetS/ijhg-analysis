# Manuscript element → script → output

| Manuscript element | Script | Output |
|---|---|---|
| Analysis dataset (82 areas × 7 waves) | `R/01_prepare_data.R` | `data/derived/analysis_dataset.rds` / `.csv` |
| Table 1 — incidence determinants | `R/04` → `R/05` | `output/tables/Table1_incidence.xlsx` |
| Table S1 — mortality | `R/04` → `R/05` | `output/tables/TableS1_mortality.xlsx` |
| Table S2 — case fatality | `R/04` → `R/05` | `output/tables/TableS2_case_fatality.xlsx` |
| Table 2 — covariates by outcome *(new)* | `config.R` (`COV_*`) | defined in code; see `docs/DATA_DICTIONARY.md` |
| Table 3 — wave definitions *(new)* | `R/02_waves_descriptives.R` | `output/tables/table3_waves.csv` |
| Table S3 — null-interval sensitivity *(new)* | `R/06_null_interval.R` | `output/tables/TableS3_null_interval.csv` |
| Table S4 — incidence with comorbidities *(new)* | `R/08_sensitivity.R` | `output/tables/TableS4_incidence_comorbid.csv` |
| Figure 1 — epidemic curve | **not automated** | legacy: `COVID_Cuba_Paper/src/grafico incidencia wave.R` |
| Figure 2 — descriptive maps | **not automated** | legacy mapping chunks |
| Figure 3 — adjusted RR maps | `R/04` → `R/09_figures.R` | `output/figures/fig3_RR_<outcome>.pdf` |
| Figure 4 — hotspot persistence | `R/04` → `R/07` → `R/09` | `output/figures/fig4_hotspots.pdf` |
| Moran's I (Results, spatial clustering) | `R/02_waves_descriptives.R` | `output/tables/morans_I.csv` |
| MCMC diagnostics *(new)* | `R/04` → `R/09` | `diagnostics.csv`, `figS_trace_acf.pdf`, `figS_diagnostics.pdf` |

## Notes

- **Mapped relative risk.** `04_fit_models.R` computes the combined adjusted RR per area
  and wave, `exp(Xβ + φ)`, expressed relative to the population-weighted city mean rate so
  that RR = 1 means "at the city average". This is the quantity mapped in Figure 3 and used
  for hotspot classification. It is *not* `exp(φ)` alone.
- **Bold / dagger conventions.** Bold-style flag = 95% CrI excludes 1.0. Dagger = interval
  ratio (upper/lower) > 20, i.e. dominated by sparse data. Both are set in `R/05_tables.R`
  (`WIDE_RATIO`).
- **Two tables are described but not produced by any script**: Figure 1 and Figure 2. If
  the revision needs them regenerated from source, they should be ported into
  `09_figures.R`; the underlying objects (`analysis_dataset.rds`) are the same.
