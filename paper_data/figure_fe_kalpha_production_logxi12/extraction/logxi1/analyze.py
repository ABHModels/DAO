from pathlib import Path
import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
RUN = ROOT / 'results/9e14af6b'
atomic_file = Path('/Users/xe26734/Downloads/cloudy/data/mewe_fluor.dat')
atomic_rows = [list(map(float, s.split()[:7])) for s in atomic_file.read_text().splitlines()
               if s.strip() and not s.startswith('#')]
atomic = np.asarray(atomic_rows)

old = np.loadtxt(RUN/'Iron/ion_fractions.dat', ndmin=2)
new = np.loadtxt(OUT/'iron_fractions_replayed.dat', ndmin=2)
assert old.shape == (48, 33) and new.shape == (48, 32)
ion_error = float(np.max(np.abs(old[:, 6:] - new[:, 5:])))
assert ion_error < 1e-5, f'Atomic replay differs from accepted ions: {ion_error}'
assert np.allclose(old[:, 3], new[:, 2], rtol=1e-7)
assert np.allclose(old[:, 4], new[:, 3], rtol=1e-5)
budget = np.loadtxt(RUN/'thermal_budget_iter63.dat')
tau, dz = budget[:, 1], budget[:, 8]
fl = np.concatenate([np.loadtxt(OUT/f'cells/fluorescence_{d}.dat', ndmin=2) for d in range(48)])
bb_files = [OUT/f'cells/bound_bound_{d}.dat' for d in range(48)]
bb = np.concatenate([np.loadtxt(p, ndmin=2) for p in bb_files
                     if any(s.strip() and not s.startswith('#') for s in p.read_text().splitlines())])
survival = np.loadtxt(OUT/'bound_bound_survival.dat', ndmin=2)
assert np.isfinite(fl).all() and np.isfinite(bb).all() and np.isfinite(survival).all()
assert np.all(fl[:, 12:] >= 0) and np.all(bb[:, 5] >= 0)
ids = fl[:, 6].astype(int)
db = atomic[ids]
assert np.all(db[:, 0] == 26)
assert np.array_equal(db[:, 1], fl[:, 7])
# DAO uses 13.6058 eV/Ryd; Cloudy's database conversion differs by 8 ppm.
assert np.allclose(db[:, 5], fl[:, 10], rtol=1e-5)
assert np.allclose(db[:, 6], fl[:, 14], rtol=1e-6)
# Fluorescence line IDs 1 and 2 are the K-alpha doublet; 3 and 4 are K-beta.
ka = (fl[:, 9] == 0) & np.isin(db[:, 4], [1, 2])
assert np.all((fl[ka, 10] >= 6200) & (fl[ka, 10] < 7000))
fka = fl[ka]
local_fl = np.zeros((48, 27))
local_grid = np.zeros((48, 27))
for r in fka:
    d, parent = int(r[0]), int(r[7])-1
    local_fl[d, parent] += r[16]
    local_grid[d, parent] += r[17]
# Bound-bound lines in the Fe K-alpha complex, including He/H-like emission.
bka = bb[(bb[:, 3] >= 6200) & (bb[:, 3] < 7000)]
lookup = {(int(r[0]), int(r[3])): r for r in survival}
local_bb = np.zeros((48, 27))
local_surv = np.zeros((48, 27))
for r in bka:
    d, ip, stage = map(int, r[:3])
    s = lookup[d, ip]
    assert np.isclose(s[5], r[5], rtol=1e-12, atol=0)
    assert 0 <= s[6] <= 1+1e-12
    local_bb[d, stage-1] += r[5]
    local_surv[d, stage-1] += s[7]

low = local_fl[:, :5].sum(axis=1)
fl_all = local_fl.sum(axis=1)
bb_all = local_bb.sum(axis=1)
total = fl_all + bb_all
injected = local_grid.sum(axis=1) + local_surv.sum(axis=1)
power_by_parent = (local_fl*dz[:, None]).sum(axis=0)
bb_by_stage = (local_bb*dz[:, None]).sum(axis=0)
low_power = float(low @ dz)
total_power = float(total @ dz)
fl_power = float(fl_all @ dz)
depth_values = np.column_stack([np.arange(48), tau, dz, new[:, 2], new[:, 3],
                               local_fl[:, :5], low, fl_all, bb_all, total,
                               low/total, low*dz, total*dz,
                               np.cumsum(low*dz)/low_power, injected])
np.savetxt(OUT/'depth_power.dat', depth_values, fmt='%.16e',
           header='depth tau_ref dz_cm T_K ne_cm-3 eps_fl_parent_FeI eps_fl_parent_FeII eps_fl_parent_FeIII eps_fl_parent_FeIV eps_fl_parent_FeV eps_low_fl eps_all_fl eps_all_bb eps_all_Kalpha low_fraction local_low_power_per_area local_total_power_per_area cumulative_low_fraction eps_material_Kalpha_after_survival\n'
                  'eps: erg cm^-3 s^-1; local_power_per_area: erg cm^-2 s^-1, all emission directions; not emergent flux. Fluorescence grouped by pre-photoionization parent; bound-bound by emitting ion.')
stages = np.arange(1, 28)
np.savetxt(OUT/'ion_power.dat', np.column_stack([stages, power_by_parent, bb_by_stage,
                                              power_by_parent/fl_power, power_by_parent/total_power]),
           fmt='%.16e', header='stage fl_Kalpha_production_per_area bb_Kalpha_production_per_area fl_fraction_of_fluorescence fl_fraction_of_all_Kalpha\n'
           'stage is photoionization parent for fluorescence, emitting ion for bound-bound. Powers are sums eps*dz, all directions, before transport.')
np.savetxt(OUT/'fluorescence_records.dat', np.column_stack([fl, db[:, 4]]), fmt='%.16e',
           header=(OUT/'cells/fluorescence_0.dat').read_text().splitlines()[0].removeprefix('# ')+' database_line_type')

peak_d = int(np.argmax(low))
summary = {
    'source_run': '9e14af6b', 'log_xi': 1, 'accepted_iteration': 63,
    'mean_intensity_from_iteration': 62,
    'definition': 'local K-alpha photon production before slab transfer, all directions; fluorescence indexed by absorbing parent ion',
    'temperature_or_transfer_rerun': False,
    'local_source_units': 'erg cm^-3 s^-1',
    'depth_integrated_source_units': 'erg cm^-2 s^-1 (generated power per slab area, not emergent flux)',
    'max_absolute_ion_fraction_replay_error': ion_error,
    'FeI_V_Kalpha_production_per_area': low_power,
    'all_Fe_Kalpha_production_per_area': total_power,
    'all_Fe_fluorescence_Kalpha_production_per_area': fl_power,
    'all_Fe_bound_bound_Kalpha_production_per_area': float(bb_all @ dz),
    'low_fraction_of_all_Kalpha': low_power/total_power,
    'low_fraction_of_fluorescence_Kalpha': low_power/fl_power,
    'FeI_V_fluorescence_power_by_parent': power_by_parent[:5].tolist(),
    'FeI_V_bound_bound_power_by_emitter': bb_by_stage[:5].tolist(),
    'fluorescence_parent_groups_fraction_of_all_Kalpha': {
        f'Fe_{lo}_{hi}': float(power_by_parent[lo-1:hi].sum()/total_power)
        for lo,hi in [(1,5),(6,16),(17,22),(23,27)]},
    'fraction_low_produced_tau_gt_1': float((low[tau>1]*dz[tau>1]).sum()/low_power),
    'fraction_low_produced_tau_gt_3': float((low[tau>3]*dz[tau>3]).sum()/low_power),
    'max_local_low_fraction': float(np.max(low/total)),
    'peak_low_emissivity_depth': peak_d,
    'peak_low_emissivity_tau': float(tau[peak_d]),
    'injected_Kalpha_source_after_survival_per_area': float(injected @ dz),
    'low_grid_fraction_of_injected_Kalpha': float((local_grid[:, :5].sum(axis=1) @ dz)/(injected @ dz)),
    'fluorescence_grid_vs_true_energy_relative_difference': float((local_grid.sum(axis=1)@dz)/fl_power-1),
}
protected = json.loads((OUT/'protected_files_before.json').read_text())
changed = [p for p,sha in protected.items() if hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=sha]
assert not changed, changed
summary['protected_files_verified'] = len(protected)
(OUT/'analysis.json').write_text(json.dumps(summary, indent=2)+'\n')
print(json.dumps(summary, indent=2))

provenance = {
    'inputs': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in
               [RUN/'params.json', RUN/'thermal_budget_iter63.dat', RUN/'moments_iter062.dat',
                RUN/'Iron/ion_fractions.dat',atomic_file]},
    'production_interface_sha256': hashlib.sha256((ROOT/'source/cloudy_interface_v2.cpp').read_bytes()).hexdigest(),
    'executable_sha256': hashlib.sha256((OUT/'diagnostic').read_bytes()).hexdigest(),
}
(OUT/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
