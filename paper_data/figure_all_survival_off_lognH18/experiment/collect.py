"""Validate and curate only a fully converged all-line P=1 atmosphere."""
from pathlib import Path
import hashlib
import json
import re
import shutil
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUN = ROOT / 'results/fe48cbde'
OUT = ROOT / 'paper_data/figure_all_survival_off_lognH18'
status = json.loads((RUN / 'thermal_status.json').read_text())
assert status['state'] == 'converged', 'Do not curate an unconverged spectrum as final'
launch = json.loads((HERE / 'launch.json').read_text())
assert launch['state'] == 'converged' and launch['production_files_unchanged'] and launch['baseline_files_unchanged']
n = status['iteration']
history = np.loadtxt(RUN / 'thermal_history.dat')
assert np.all(np.isfinite(history)) and history[-1, 0] == n
last = history[-3:]
assert len(last) == 3 and np.all(last[:, 1:4] < .003)
assert np.all(last[:, 4] < .01) and np.all(np.abs(last[:, 5] - 1) < .01)
assert np.all(np.abs(history[:, 6]) <= 1e-4)

lines = np.loadtxt(RUN / f'line_escape_lines_iter{n:03d}.dat', usecols=range(15))
emitting = lines[:, 4] > 0
assert emitting.any() and np.allclose(lines[emitting, 7], 1, rtol=0, atol=1e-9)
assert np.all(lines[:, 6] == 0) and np.all(lines[:, 8] == 0)

for name in ('input','validation','experiment','output/pdf','data'):
    (OUT / name).mkdir(parents=True, exist_ok=True)


def upper_flux(path):
    a = np.loadtxt(path)
    mu, wt = np.polynomial.legendre.leggauss(8)
    header = path.read_text().splitlines()[:12]
    recorded_mu = [float(re.search(r'mu=([^)]*)', line).group(1)) for line in header if 'mu=' in line]
    assert np.allclose(mu, recorded_mu, atol=5e-7, rtol=0)
    positive = mu > 0
    fe = 2 * np.pi * (a[:, 3:][:, positive] @ (wt[positive] * mu[positive])) * 1e3
    energy = a[:, 0] / 1e3
    assert np.all(np.isfinite(fe)) and np.all(fe >= 0)
    return np.column_stack((energy, fe, energy * fe)), a[:, 1:3]


std, standard_incident = upper_flux(ROOT / 'results/5b7fe321/emergent_iter019.dat')
new, new_incident = upper_flux(RUN / f'emergent_iter{n:03d}.dat')
assert np.array_equal(std[:, 0], new[:, 0])
assert np.array_equal(standard_incident, new_incident)
for name, data in (('DAO_standard', std), ('DAO_all_P1', new)):
    np.savetxt(OUT / 'input' / f'{name}.dat', data, fmt='%.16e',
               header='E_keV F_E_erg_cm-2_s-1_keV-1 E_F_E_erg_cm-2_s-1; upper outgoing flux = 2*pi*sum(w*mu*I), mu>0')

for name in ('reflionx_flux_cgs.dat','xillverCp_angle_average.dat'):
    shutil.copy2(ROOT / 'paper_data/figure_oxygen_highdensity/input' / name, OUT / 'input' / name)
for src, dest in ((RUN / f'thermal_budget_iter{n}.dat', 'thermal_budget_all_P1.dat'),
                  (ROOT / 'results/5b7fe321/thermal_budget_iter19.dat', 'thermal_budget_standard.dat'),
                  (RUN / f'profile_iter{n:03d}.dat', 'profile_all_P1.dat'),
                  (ROOT / 'results/5b7fe321/profile_iter019.dat', 'profile_standard.dat'),
                  (RUN / 'thermal_history.dat', 'thermal_history.dat'),
                  (RUN / 'params.json', 'params.json'),
                  (RUN / 'O/ion_fractions.dat', 'oxygen_ions_all_P1.dat'),
                  (RUN / 'Iron/ion_fractions.dat', 'iron_ions_all_P1.dat')):
    shutil.copy2(src, OUT / 'input' / dest)
for name in ('launch.json','build_manifest.json','verification.log'):
    shutil.copy2(HERE / name, OUT / 'validation' / name)
for name in ('prepare.py','experiment.mk','experiment.patch','run.py','collect.py'):
    shutil.copy2(HERE / name, OUT / 'experiment' / name)
shutil.copytree(HERE / 'source', OUT / 'experiment/source', dirs_exist_ok=True)
shutil.copy2(HERE / 'plot_comparison.py', OUT / 'plot_comparison.py')
shutil.copy2(HERE / 'README.md', OUT / 'experiment/README.md')

peak_region = np.flatnonzero((std[:,0] >= .64) & (std[:,0] <= .67))
old_peak = peak_region[np.argmax(std[peak_region,2])]
new_peak = peak_region[np.argmax(new[peak_region,2])]
old_T = np.loadtxt(OUT / 'input/thermal_budget_standard.dat')[:,2]
new_T = np.loadtxt(OUT / 'input/thermal_budget_all_P1.dat')[:,2]
report = dict(
    baseline='5b7fe321', run='fe48cbde', final_iteration=n,
    scope='All bound-bound survival factors P=1; full temperature, ion populations, opacities and radiation recomputed from normal initialization.',
    wall_seconds=launch['wall_seconds'],
    all_emitting_lines_have_P1=True, emitting_line_count=int(emitting.sum()),
    max_dlogT_dex=float(last[-1,1]), max_relative_dJ=float(last[-1,2]),
    max_local_residual=float(last[-1,3]), abs_volume_residual_over_Fin=float(last[-1,4]),
    Fout_over_Fin=float(last[-1,5]), boundary_energy_error_percent=float(100*abs(last[-1,5]-1)),
    max_identity_error=float(np.max(np.abs(history[:,6]))),
    last_three_convergence_rows=last.tolist(),
    oviii_full_spectrum_peak=dict(standard_energy_keV=float(std[old_peak,0]),
        all_P1_energy_keV=float(new[new_peak,0]), standard_absolute_EF_E=float(std[old_peak,2]),
        all_P1_absolute_EF_E=float(new[new_peak,2]),
        absolute_peak_change_percent=float(100*(new[new_peak,2]/std[old_peak,2]-1))),
    temperature_K=dict(standard_surface=float(old_T[0]), all_P1_surface=float(new_T[0]),
                       standard_bottom=float(old_T[-1]), all_P1_bottom=float(new_T[-1])),
    input_sha256={str(p.relative_to(OUT)):hashlib.sha256(p.read_bytes()).hexdigest()
                  for p in sorted((OUT/'input').glob('*'))},
)
(OUT / 'measurements.json').write_text(json.dumps(report,indent=2)+'\n')
(OUT / 'README.md').write_text('''# Fully re-equilibrated all-line P=1 comparison

The baseline is DAO run 5b7fe321. The comparison is run fe48cbde with the
additional DAO survival correction disabled for ALL bound-bound lines.
Temperature, ion populations, opacities, emissivities and radiation are solved
again from normal initialization with unchanged convergence criteria.
This is distinct from the fixed-atmosphere O VIII-only test.

The incident parameters are log xi=3, log nH=18, Gamma=1.8, kTe=60 keV,
kTbb=0.01 keV, reflionx flux normalization, solar Fe, 48 depth cells,
reference Thomson depth 5, angle-averaged Compton scattering.

Run `python plot_comparison.py` here to recreate the figure using only the
curated inputs (NumPy and Matplotlib required). All panels show EF_E and retain
native energy sampling without smoothing. BOTH DAO curves share the STANDARD
DAO normalization in each panel. XillverCp and reflionx retain their individual
normalizations from the previous comparison. Absolute outgoing fluxes are
stored in input/DAO_standard.dat and input/DAO_all_P1.dat, before normalization.
The units and angular-flux definition are recorded in their headers.

The line probabilities in this artificial P=1 experiment are imposed rather
than inferred from its optical depths. Continuum absorption and Compton
scattering remain active. Collisional rates in Cloudy's statistical equilibrium
are retained. Fluorescence and the existing H I Lyalpha thin-reference
normalization remain unchanged.

measurements.json records the final convergence criteria, total runtime,
absolute O VIII peak change and input hashes. Changes in the full line peak
include changes in the atmosphere and other emission; they are not an isolated
measurement of O VIII photon survival at fixed structure.

The experiment directory records the exact source patch, isolated build and
launch scripts; validation contains source/input provenance and the two-path
line-source verification. Repeating the thermal solve requires the original
DAO/Cloudy/HEASoft installation, kernel and source versions. Production source
files and the original accepted results were verified unchanged.
''')
print(json.dumps(report,indent=2))
