# Area-level factors associated with COVID-19 in Havana, Cuba

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23857423.svg)](https://doi.org/10.5281/zenodo.23857423)

Analysis code for the manuscript *Area-level factors associated with COVID-19 incidence,
mortality and case fatality across epidemic waves in Havana, Cuba: a wave-specific Bayesian
spatial ecological study*.

The study covers the 82 health areas of Havana between 11 March 2020 and 31 May 2022,
divided into seven epidemic waves. Three outcomes are modelled separately for each wave —
cumulative incidence, cumulative mortality and case fatality — with Bayesian Poisson
regression and a Leroux conditional autoregressive prior for the spatial random effect.

**Seventeen models** are fitted:

| Outcome | Waves | Models |
|---|---|---|
| Incidence | Early, Wuhan, D614G, Beta, Delta, Omicron BA.1, Omicron BA.2 | 7 |
| Mortality | Early, Wuhan, D614G, Beta, Delta | 5 |
| Case fatality | Early, Wuhan, D614G, Beta, Delta | 5 |

The two Omicron waves are excluded from the severity models because they recorded 4 and 2
deaths respectively in the whole city. Early (48 deaths) and Wuhan (20 deaths) are below
`MIN_DEATHS_SEVERITY` and are fitted but flagged `sparse_outcome`: with 18 covariates they
are over-parameterised and their estimates are exploratory.

## Requirements

R ≥ 4.1 with: CARBayes (6.x), spdep, sf, coda, dplyr, tidyr, tibble, readxl, writexl,
scales, ggplot2, gridExtra, MASS. The analyses reported in the paper were run with R 4.1.3
and CARBayes 6.1.1 on the lab workstation.

## How to run

```bash
Rscript run_all.R                      # main analysis (17 models, tables, figures)
Rscript run_all.R --with-sensitivity   # plus the reviewer-requested analyses
QUICK_RUN=1 Rscript run_all.R          # short chains, smoke test only
```

`run_all.R` writes a manifest of every output and a copy of `sessionInfo()` to
`output/logs/`.

Alternative runs keep their outputs in a sub-folder via `RUN_TAG`:

```bash
HACIN_SOURCE=fichero RUN_TAG=hacin_fichero Rscript run_all.R   # alternative crowding definition
```

### Running the 17 models across processes

`tools/run_parallel.R` spreads the model fits over independent processes and writes the same
tables. The sensitivity analyses have their own driver:

```bash
N_WORKERS=8  Rscript tools/run_parallel.R
N_WORKERS=8  Rscript tools/run_sens_parallel.R
SECTIONS=C N_WORKERS=8 Rscript tools/run_sens_parallel.R   # one section only
```

## Structure

```
config.R                  paths, constants, covariate sets, prior/MCMC settings
run_all.R                 driver
R/
  00_setup.R              packages and helpers
  01_prepare_data.R       corrected panel + shapefiles -> analysis dataset (82 x 7)
  02_waves_descriptives.R wave periods; Moran's I
  03_standardise.R        covariate scaling
  04_fit_models.R         the 17 wave-specific Leroux CAR models
  05_tables.R             Tables 1, S1, S2
  06_null_interval.R      Table S3 (criterion sensitivity)
  07_hotspots.R           hotspot classification with and without posterior uncertainty
  08_sensitivity.R        comorbidities, binomial CFR, weights, priors, collinearity, waves
  09_figures.R            Figures 3 and 4; MCMC diagnostics
  14_vaccination_severity.R   vaccination in the Delta severity models
  16_figures_descriptive.R    Figure 2 and the continuous covariates
  17_figure1.R            Figure 1 (epidemic curve and interventions)
  18_regularised_severity.R   regularised severity models
  19_figure_S1.R          Figure S1 (covariate maps)
  20_severity_incidence_adjustment.R  severity adjusted for same-wave incidence
  21_collinearity.R       correlations and VIF
  22_wave_boundaries.R    wave-boundary sensitivity
  23_table2_covariates.R  Table 2 (covariates by model)
  24_negative_binomial.R  negative-binomial sensitivity (non-spatial)
data/raw/                 restricted inputs (not included; see data/README.md)
data/synthetic/           synthetic test dataset (included)
docs/                     corrections, data dictionary, script->output map, reviewer crosswalk
output/tables/            results
output/figures/           figures
```

## Model specification

For health area *i* in wave *j*:

```
log(mu_ij) = log(E_ij) + alpha + x_ij'beta + phi_i
```

* `E_ij` is the offset: population for incidence and mortality, confirmed cases for case
  fatality.
* `phi_i` is a spatially structured random effect with a Leroux prior, as implemented in
  `S.CARleroux` (CARBayes):

```
phi_i | phi_-i ~ Normal( rho * sum_j w_ij phi_j / (rho * sum_j w_ij + 1 - rho),
                         tau2 / (rho * sum_j w_ij + 1 - rho) )
```

with `w_ij` from first-order queen contiguity, `rho` the spatial dependence parameter and
`tau2` the variance parameter.

**Covariates.** Incidence models carry 12 structural covariates: population density,
education, ageing index, overcrowding, shanty towns, distance to the nearest family medical
centre, population served per centre, CMF density, agglomeration sites, rental housing
density, piped-water access and precarious housing. Severity models drop density (it takes
the VIF above 5 once the comorbidities enter) and add six chronic-disease prevalence
variables; mortality is additionally adjusted for the same-wave incidence. Vaccination
coverage is included in the Delta incidence model only.

All covariates are rescaled to [0, 1] across the 82 areas, so every reported relative risk
refers to the full observed range of the covariate.

**Estimation.** Three chains per model, 1M–6M iterations for incidence and 150k for
severity, with convergence assessed by the Gelman–Rubin statistic and effective sample size.
Primary models use the CARBayes default priors (beta ~ Normal(0, 100000); precision of the
CAR effect ~ Gamma(1, 0.01), i.e. inverse-gamma(1, 0.01) on tau2; rho uniform on [0, 1]).
The priors declared in the submitted Methods, and heavier and lighter alternatives, are
fitted as a sensitivity analysis in `08_sensitivity.R`.

## Data and synthetic test dataset

The health-area-level surveillance data are held by the Cuban Ministry of Public Health and
are **not** included in this repository. The analysis reads a corrected panel in which
several columns of the original file were re-aligned with the health-area codes; the
corrections and the audit are documented in `docs/CORRECTIONS.md`, and `data/README.md`
lists the files the pipeline expects.

To make the pipeline runnable without the restricted data, `data/synthetic/` contains a
**synthetic test dataset** with the same structure (82 areas x 7 waves, same files and
columns), generated by `tools/make_synthetic_data.R`. The geography is a Voronoi
tessellation of a 30 x 24 km rectangle and every covariate and outcome is simulated, so
nothing in it corresponds to a real health area or a real person. `config.R` uses the real
data when the corrected panel is present in `data/raw/` and falls back to the synthetic
dataset otherwise, printing a message when it does. **Results produced in that mode are
illustrative only and must not be interpreted.**

```bash
Rscript tools/make_synthetic_data.R   # regenerate the synthetic inputs
```

## Citation

If you use this code, please cite the article and the archived release. `CITATION.cff` also
lets GitHub display a "Cite this repository" button.

The archived version is **v1.0.0**, with concept DOI
[10.5281/zenodo.23857423](https://doi.org/10.5281/zenodo.23857423) (this version:
[10.5281/zenodo.23857424](https://doi.org/10.5281/zenodo.23857424)).

> Sanchez-Valdes L, Sauleau EA, Remond Noa R, Palomino Hernandez JL, Veranes Miranda A,
> Handschumacher P, Van der Stuyft P, Mas Bermejo P. Area-level factors associated with
> COVID-19 incidence, mortality and case fatality across epidemic waves in Havana, Cuba: a
> wave-specific Bayesian spatial ecological study. *International Journal of Health
> Geographics*. [volume, year, DOI to be added on publication]

## License

The code is released under the MIT License; see `LICENSE`. The health-area surveillance
data are not included in this repository and are not covered by this license; they remain
under the custodianship of the Cuban Ministry of Public Health.

## Contact

Lizet Sanchez-Valdes, Clinical Research Division, Center of Molecular Immunology, Havana,
Cuba — lsanchez@cim.sld.cu
