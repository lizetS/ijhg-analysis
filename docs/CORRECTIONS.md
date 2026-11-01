# Corrections to the input data (supersedes earlier notes on population and covariates)

The original panel `habana_covid_integration_by_variant_Wave_310825.xlsx` has several
covariate columns attached to the wrong health area. The analysis now reads
`data/raw/habana_panel_corregido.xlsx` (sheets: `panel_corregido`, `cambios`, `pendientes`).
Audit: `python3 tools/audit_panel_alignment.py`.

## What was wrong in the original panel
- `poblacion_V`, `DensPob`, `DensPob_C`: not aligned with the area (log-correlation with cases
  0.09 vs 0.36 for `poblacion`, first wave). 57 of 82 areas differ from `poblacion`.
- The manuscript's "persistent hotspots" (Figure 4) are reproduced exactly (82/82 areas) by counting
  the waves with `RR_I_V` >= 1.5. `RR_I_V` is the crude incidence ratio computed with `poblacion_V`.
  All 11 persistent areas have a `poblacion_V` of 33-69 % of their real population, so they appear
  inflated in every wave. With a correct denominator no area is persistent.
- `Prop_Hip/Diab/Card/Epoc/Asma` and `Dist_x_AS` were a row permutation of the values in
  `AS_info_ok_todas_var.xlsx` (62 of 82 areas carried another area's values), and that file is in turn
  misaligned with `informacion_comorbilidades.xlsx` (19 of 81 areas in the right place).

## What the corrected panel uses
| Variable | Source / rule |
|---|---|
| `poblacion` | panel column `poblacion` (decided by the author) |
| density | `poblacion` / area from the geometry (`density_corrected` in the panel) |
| `Prop_*` (6 comorbidities) | counts in `informacion_comorbilidades.xlsx` / `poblacion` x 100, by area |
| `Dist_x_AS`, `Dis_x_CMF_A`, `ageing_index_cub`, `mayores60` | `AS_info_ok_todas_var.xlsx`, by Etiqueta |
| `hacin` | panel/census definition (single definition; the file-derived column was removed by decision of the author) |
| `vacunados` | Delta wave only |
| removed | `poblacion_V`, `DensPob`, `DensPob_C`, all `RR_*`, `fitted_*`, `hotspot_type`, `Inc` |

## Open points (also in the `pendientes` sheet)
1. Crowding: two definitions (census estimate vs health-area report); which one the manuscript
   used is unknown. Run both: `HACIN_SOURCE=panel RUN_TAG=hacin_panel Rscript run_all.R`.
2. `Dens_x_CMF` (panel) vs `Dens_x_CMF_A` (reviewed file, kept as `Dens_x_CMF_fichero`): the analysis uses `Dens_x_CMF`.
3. Nine areas where `poblacion` differs between the panel and the reviewed file: Calabazar, Capdevila,
   Santiago de las Vegas, Wajay, Efrain Mayor, Rafael Valdes, Gregorio Valdes, Rafael Zilueta, 14 de Junio.
4. `Dis_x_CMF_A` of Betancourt Neninger: corrected to 1726.13 (decision of the author).
5. `informacion_comorbilidades.xlsx`: Capdevila 60+ = 1,408 (reviewed file: 7,537); Mario Escalona and
   13 de Marzo both 7,012. Asthma maximum 33.9 % in one area.
6. `GradE`, `PcAAgua`, `PcvivD` come from the auxiliary shapefile and were not checked against the reviewed file.
7. All model results produced before this correction (including `output/tables/*.csv` and the
   statements in `REPRODUCTION_FINDINGS.md`, `REVISION_TEXT.md`, `NOTE_TO_COAUTHORS.md`) used the
   uncorrected inputs and must be regenerated.

## Severity waves
Mortality and case-fatality models are now fitted for Early, Wuhan, D614G, Beta and Delta
(`SEVERITY_WAVES`). Omicron BA.1 and BA.2 (4 and 2 deaths) are excluded. Early (48 deaths) and Wuhan
(20 deaths) are below `MIN_DEATHS_SEVERITY` and are flagged `sparse_outcome`: with 18 covariates they are
overparameterised, so their wide intervals should be reported and the results presented as exploratory.
Text in `REVISION_TEXT.md` and `Response_to_Reviewers_v2.md` that says "D614G, Beta and Delta" for severity is obsolete.

## Priors: declared versus used
- Submitted manuscript (Methods): "weakly informative normal priors" on the covariate effects and a "gamma
  prior" on the precision of the CAR component. No numbers are given; rho is not mentioned.
- Earlier README / draft Methods: N(0, 100) for beta and Gamma(0.5, 0.0005) for tau2.
- Code (`fit_model` passes no priors for the primary models), so CARBayes defaults apply (CARBayes
  documentation): beta ~ N(0, variance 100,000); tau2 ~ Inverse-Gamma(shape 1, scale 0.01), which is
  Gamma(1, rate 0.01) on the precision; rho estimated with a uniform prior on [0, 1].
- Therefore the primary models do NOT use the priors of the README / draft Methods. The family
  (normal, gamma on the precision) matches the manuscript; the parameter values do not, and
  N(0, 1e5) is vague rather than "weakly informative" on the min-max scale.
- Sensitivity: `08_sensitivity.R` section D now also fits the declared README priors
  (`declared_N100_G05`) and N(0, 100) on beta alone. `18_regularised_severity.R` adds N(0, 1) on the z scale.
- Software: the results in `archive_uncorrected` were produced with R 4.6.0 and CARBayes 6.1.1; the
  manuscript states R 4.1.1 and CARBayes 5.3. Report the versions actually used.
