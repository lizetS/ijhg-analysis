# Reviewer point → where it is answered in code

Reviewer numbering follows the files supplied to the authors: **R1** = section-by-section
review (`reviewer 1.docx`, internally headed "Reviewer 2"), **R2** = 17-point review.

| Point | Topic | Implemented in | Output |
|---|---|---|---|
| R1 1.11 | Incidence models with comorbidities | `R/08` section A | `TableS4_incidence_comorbid.csv`, `TableS4_comparison.csv` |
| R1 1.12, R2 2.4 | Covariate standardisation | `R/03` | `scaling_parameters.csv` |
| R1 1.13, 1.17, R2 2.17 | Leroux specification / software | `R/04` (`S.CARleroux`) | `diagnostics.csv` (rho, tau2) |
| R1 1.14, R2 2.6 | Chains and convergence | `R/04` (3 chains, R̂, ESS, ACF) | `diagnostics.csv`, `figS_trace_acf.pdf`, `figS_diagnostics.pdf` |
| R1 1.15 | Zero deaths / zero cases | `R/01` | `zero_event_counts_by_wave.csv` |
| R1 1.16, 1.21 | Meaning of mapped RR | `R/04` (combined adjusted RR), `R/09` | `area_rr_summary.csv`, `fig3_RR_*.pdf` |
| R1 1.18 | Null-interval sensitivity | `R/06` | `TableS3_null_interval.csv` |
| R1 1.19 | Hotspot counts and percentages | `R/07` | `hotspot_counts.csv` |
| R1 1.23 | Multiple comparisons | `R/04` (`pr_gt_1`), `R/06` | `coefficients.csv`, `TableS3_summary.csv` |
| R2 2.3, 2.7 | Binomial case-fatality model | `R/08` section B | `sensitivity_binomial_cfr.csv` |
| R2 2.5 | Collinearity | `R/08` section E | `collinearity.csv` |
| R2 2.8 | CAR variance prior | `R/08` section D | `sensitivity_priors.csv` |
| R2 2.9 | Spatial weights | `R/08` section C | `sensitivity_weights.csv` |
| R2 2.10 | Hotspot uncertainty | `R/07` | `hotspot_classification.csv`, `hotspot_reclassified.csv` |
| R2 2.14 | Vaccination measure | `R/01` (`cobertura_vac`), `config.R` | `scaling_parameters.csv` (n = 82) |
| R2 2.15 | Wave definition robustness | `R/08` section F | `wave_robustness.csv` |
| R2 2.16 | Very large / imprecise estimates | `R/05` (dagger flag), `R/01` | Tables S1/S2, `zero_event_counts_by_wave.csv` |

## Still to be added by hand

- **99% credible intervals** for the interpreted associations (R1 1.23). The posterior
  object in `04` keeps only a thinned trace; to report 99% CrI, extend `fit_model()` to
  compute `quantile(pooled, c(0.005, 0.995))` alongside the 95% interval.
- **Delta mortality and case-fatality models including vaccination** (optional, R2 2.14).
  One call each: `fit_model("mortality", "Delta", c(covariates_for("mortality"), "cobertura_vac"))`.
- **Testing intensity / positivity** (R2 2.13) cannot be implemented: the data are not
  available at health-area level. Handled in the text, not in code.
