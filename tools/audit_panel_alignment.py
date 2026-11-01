"""Audit of the corrected panel (run from the repository root).
Checks, without R: required columns, missing values, constancy of area-level
covariates across waves, agreement of area with the shapefile geometry, and the
original-panel problems that motivated the correction.
Usage: python3 tools/audit_panel_alignment.py
"""
import numpy as np, pandas as pd, shapefile
from shapely.geometry import shape

RAW = 'data/raw/'
new = pd.read_excel(RAW + 'habana_panel_corregido.xlsx', sheet_name='panel_corregido')
old = pd.read_excel(RAW + 'habana_covid_integration_by_variant_Wave_310825.xlsx')
old = old[old.covid_variant.notna() & old.ASalud.notna()]

r = shapefile.Reader(RAW + 'shapefiles/AS_vuln_vs_casos_8_periodos')
shp = pd.DataFrame(r.records(), columns=[x[0] for x in r.fields[1:]])
shp['area_km2'] = [shape(s.__geo_interface__).area / 1e6 for s in r.shapes()]
ra = shapefile.Reader(RAW + 'shapefiles/otras variables/AS_con_enlace_lista')
aux = pd.DataFrame(ra.records(), columns=[x[0] for x in ra.fields[1:]])

ok = True
def check(name, cond, detail=''):
    global ok
    ok &= bool(cond)
    print(('PASS ' if cond else 'FAIL ') + name + (' -- ' + str(detail) if detail != '' else ''))

check('574 rows = 82 areas x 7 waves', len(new) == 574 and new.codigo_as.nunique() == 82)
check('no duplicated area x wave', not new.duplicated(['codigo_as', 'covid_variant']).any())
check('all codes in both shapefiles', set(new.codigo_as) == set(shp.Etiqueta) == set(aux.Etiqueta))

static = ['poblacion', 'areaKm2', 'ageing_index_cub', 'hacin', 'CantDBprec', 'Dis_x_CMF_A', 'PobAt_CMF_A',
          'Aglomeraci', 'casasAlq', 'Dens_x_CMF', 'Prop_Hip', 'Prop_Diab', 'Prop_Card', 'Prop_Epoc',
          'Prop_cancer', 'Prop_Asma', 'Dist_x_AS']
for c in static:
    check(f'{c}: present, complete, constant within area',
          c in new and new[c].notna().all() and (new.groupby('codigo_as')[c].nunique() <= 1).all())

w = new[new.covid_variant == 'Wuhan1'].merge(shp[['Etiqueta', 'area_km2']], left_on='codigo_as', right_on='Etiqueta')
check('areaKm2 equals shapefile geometry area', np.allclose(w.areaKm2, w.area_km2, rtol=1e-3))
check('vacunados only in Delta', new.groupby('covid_variant').vacunados.apply(lambda s: s.notna().sum()).to_dict()
      == {'Beta': 0, 'D614G': 0, 'Delta': 82, 'Omicron1': 0, 'Omicron2': 0, 'Wuhan1': 0, 'Wuhan2': 0})
for c in ['poblacion_V', 'DensPob', 'DensPob_C', 'RR_I_V', 'hotspot_type']:
    check(f'{c} removed from the corrected panel', c not in new.columns)

# cases/deaths unchanged with respect to the original panel
m = new.merge(old[['codigo_as', 'covid_variant', 'cases_per_variant', 'deaths_per_variant']],
              on=['codigo_as', 'covid_variant'], suffixes=('', '_old'))
check('cases and deaths identical to original panel',
      (m.cases_per_variant == m.cases_per_variant_old).all() and (m.deaths_per_variant == m.deaths_per_variant_old).all())

# the problems found in the ORIGINAL panel (informational)
q = old[old.covid_variant == 'Wuhan1']
print('\nOriginal panel, informational:')
print(' corr(log poblacion_V, log cases) = %.2f vs corr(log poblacion, log cases) = %.2f' % (
    np.corrcoef(np.log(q.poblacion_V), np.log(q.cases_per_variant))[0, 1],
    np.corrcoef(np.log(q.poblacion), np.log(q.cases_per_variant))[0, 1]))
n15 = (old.RR_I_V >= 1.5).groupby(old.codigo_as).sum()
print(' areas classed persistent hotspot (RR_I_V >= 1.5 in >= 5 waves): %d' % (n15 >= 5).sum())
print('\nOVERALL:', 'ALL CHECKS PASSED' if ok else 'SOME CHECKS FAILED')
