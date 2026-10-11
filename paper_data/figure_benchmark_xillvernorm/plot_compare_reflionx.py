"""Figure 3 — DAO vs reflionx vs xillvercp over a log xi scan.

Author:      Yimin Huang
Affiliation: Fudan University; University of Bristol
Email:       huangym23@m.fudan.edu.cn
"""

import argparse
import hashlib
import re
import json
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'dao-benchmark-matplotlib'))
import numpy as np
import matplotlib as mpl
if '--show' not in sys.argv:
    mpl.use('Agg')
import matplotlib.pyplot as plt

# ── Nature-style rcParams ────────────────────────────────────────────────────
mpl.rcParams.update({
    'font.family':       'sans-serif',
    'font.sans-serif':   ['Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size':          7,
    'axes.labelsize':     8,
    'axes.titlesize':     8,
    'axes.linewidth':     0.6,
    'xtick.labelsize':    7,
    'ytick.labelsize':    7,
    'xtick.direction':    'in',
    'ytick.direction':    'in',
    'xtick.top':          True,
    'ytick.right':        True,
    'xtick.major.size':   3.0,
    'ytick.major.size':   3.0,
    'xtick.minor.size':   1.8,
    'ytick.minor.size':   1.8,
    'xtick.major.width':  0.6,
    'ytick.major.width':  0.6,
    'xtick.minor.width':  0.5,
    'ytick.minor.width':  0.5,
    'legend.fontsize':    6.5,
    'legend.frameon':     False,
    'legend.handlelength': 1.6,
    'legend.handletextpad': 0.5,
    'pdf.fonttype':       42,
    'ps.fonttype':        42,
})

# ── Upper-face outward normal flux ─────────────────────────────────────────
def load_emergent_flux(em_file):
    """Return E [eV], 2*pi*sum(w*mu*I) [erg/cm^2/s/eV], and |mu_inc|."""
    mu_v = []
    with open(em_file) as f:
        for line in f:
            if not line.startswith('#'):
                break
            m = re.search(r'I\(mu=([-\d.]+)\)', line)
            if m:
                mu_v.append(float(m.group(1)))
    data = np.loadtxt(em_file, comments='#')
    E = data[:, 0]
    mu_a, wt_a = np.polynomial.legendre.leggauss(len(mu_v))
    assert data.shape[1] == len(mu_a) + 3
    assert np.allclose(mu_a, mu_v, rtol=0, atol=5e-7)
    assert np.all(np.isfinite(data)) and np.all(np.diff(E) > 0)
    outgoing = mu_a > 0
    F = 2 * np.pi * np.sum(data[:, 3:][:, outgoing] *
                          (wt_a[outgoing] * mu_a[outgoing]), axis=1)
    assert np.all(F >= 0)
    incoming = np.flatnonzero(mu_a < 0)
    illuminated = incoming[np.argmax(np.sum(data[:, 3:][:, incoming], axis=0))]
    return E, F, abs(mu_a[illuminated])


# Self-contained: every input file sits next to this script.
HERE = os.path.dirname(os.path.abspath(__file__))


def load_run(record):
    """Load DAO emergent spectrum, reflionx, xillver, params — all local."""
    hash_ = record['DAO_hash']
    em_file = os.path.join(HERE, record['DAO_file'])
    print(f'[{hash_}] final iteration {record["DAO_iteration"]}: {os.path.basename(em_file)}')
    E_eV, F_up, mu_inc = load_emergent_flux(em_file)

    with open(os.path.join(HERE, record['DAO_params'])) as f:
        par = json.load(f)

    reference = record['reference_parameters']
    for key in ('zeta', 'corona', 'nh', 'Gamma', 'kT_e', 'kT_bb', 'Afe'):
        assert par[key] == reference[key], (hash_, key)
    assert par['zeta'] == record['log_xi'] and record['DAO_state'] == 'converged'

    xv = np.loadtxt(os.path.join(HERE, record['xillver_file']))
    E_xv = xv[:, 0] * 1e3                    # keV -> eV
    EFE_xv = xv[:, 2] / E_xv                 # column 2 is E·F_E → divide for F_E

    rf = np.loadtxt(os.path.join(HERE, record['reflionx_file']))
    E_rf = rf[:, 0] * 1e3
    EFE_rf = rf[:, 1] / E_rf

    return dict(par=par, E=E_eV, F=F_up, mu_inc=mu_inc, record=record,
                E_xv=E_xv, F_xv=EFE_xv,
                E_rf=E_rf, F_rf=EFE_rf)


def interp_log(E_target, E_src, F_src):
    return np.exp(np.interp(np.log(E_target),
                            np.log(E_src),
                            np.log(np.clip(F_src, 1e-300, None))))


# ── Runs to plot, ordered by increasing log ξ ────────────────────────────────
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--show', action='store_true', help='Open an interactive zoom/pan window')
args = parser.parse_args()
with open(os.path.join(HERE, 'benchmark_inputs.json')) as f:
    records = json.load(f)['runs']
records = sorted(records, key=lambda r: r['log_xi'])
assert [r['log_xi'] for r in records] == [1, 2, 3]
PANEL_LABELS = ['a', 'b', 'c']
runs = [load_run(r) for r in records]
trapz = getattr(np, 'trapezoid', None) or np.trapz
plot_records = []

# Normalisation band: 20–50 keV continuum (Compton-hump, line-free).
# Each model is scaled so that ⟨F_E⟩ over this band equals 1, isolating
# spectral *shape* differences (lines/edges/soft excess) from absolute-flux
# differences between code conventions.
E_LO, E_HI = 20e3, 50e3

# Nature-friendly palette (Wong, colour-blind safe)
C_DAO   = '#0072B2'   # blue
C_REFL  = '#444444'   # dark grey
C_XILL  = '#D55E00'   # vermillion

# ── Figure ───────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.5), sharey=True)

for ax, run, label in zip(axes, runs, PANEL_LABELS):
    E   = run['E']
    F   = run['F']
    F_xv = interp_log(E, run['E_xv'], run['F_xv'])
    F_rf = interp_log(E, run['E_rf'], run['F_rf'])

    band = (E >= E_LO) & (E <= E_HI)
    dE_band = E_HI - E_LO
    # Mean F_E over the 20–50 keV continuum band
    mean_dao = trapz(F[band],    E[band]) / dE_band
    mean_xv  = trapz(F_xv[band], E[band]) / dE_band
    mean_rf  = trapz(F_rf[band], E[band]) / dE_band
    assert all(np.isfinite(v) and v > 0 for v in (mean_dao, mean_xv, mean_rf))

    F_n     = F    / mean_dao
    F_xv_n  = F_xv / mean_xv
    F_rf_n  = F_rf / mean_rf
    assert all(np.all(np.isfinite(v)) for v in (F_n, F_xv_n, F_rf_n))
    files = [run['record'][key] for key in
             ('DAO_file', 'DAO_params', 'reflionx_file', 'xillver_file')]
    plot_records.append({**run['record'], 'mu_inc_actual': float(run['mu_inc']),
                         'band_means': {'DAO': float(mean_dao), 'reflionx': float(mean_rf),
                                        'xillvercp': float(mean_xv)},
                         'input_sha256': {name: hashlib.sha256(Path(HERE, name).read_bytes()).hexdigest()
                                          for name in files}})

    # Shade the normalisation band
    ax.axvspan(E_LO, E_HI, color='0.85', alpha=0.5, lw=0, zorder=0)

    ax.loglog(E, F_n,    color=C_DAO,  lw=0.55, alpha=0.85,
              label='DAO', zorder=3)
    ax.loglog(E, F_rf_n, color=C_REFL, lw=0.9,  alpha=0.95,
              label='reflionx', zorder=5)
    ax.loglog(E, F_xv_n, color=C_XILL, lw=0.9,  alpha=0.95, ls='--',
              label='xillvercp',  zorder=5)

    ax.set_xlim(1e2, 1e6)
    # Include the new soft-energy peaks without clipping them at 10^3.
    ax.set_ylim(1e-5, 1e4)
    ax.set_xlabel(r'Energy (eV)')

    # Panel label (top-left)
    ax.text(0.03, 0.96, label, transform=ax.transAxes,
            fontsize=9, fontweight='bold', va='top', ha='left')

    # ξ annotation (top-right inside panel)
    ax.text(0.97, 0.96, rf'$\log\xi = {run["par"]["zeta"]}$',
            transform=ax.transAxes,
            fontsize=7, va='top', ha='right')

    ax.tick_params(which='both')

axes[0].set_ylabel(r'$\hat{F}_E$ (normalised)')
axes[0].legend(loc='lower left', fontsize=6.5, borderpad=0.3)

# Normalisation note in the empty lower-left region of the last panel.
axes[-1].text(
    0.03, 0.05,
    r'$\hat{F}_E = F_E \, / \, \langle F_E\rangle_{20-50\,\mathrm{keV}}$',
    transform=axes[-1].transAxes,
    ha='left', va='bottom', fontsize=6.5,
    bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='0.7', lw=0.4))

# Incident-spectrum parameters as in-figure text (panel a, bottom-left)
par0 = runs[0]['par']
incident_txt = (
    'nthcomp\n'
    rf'$\Gamma = {par0["Gamma"]}$' + '\n'
    rf'$kT_e = {par0["kT_e"]}\,$keV' + '\n'
    rf'$kT_{{\rm bb}} = {par0["kT_bb"]}\,$keV' + '\n'
    rf'$n_{{\rm H}} = 10^{{{par0["nh"]}}}\,$cm$^{{-3}}$' + '\n'
    rf'$|\mu_{{\rm inc}}| = {runs[0]["mu_inc"]:.4f}$'
)
axes[1].text(
    0.03, 0.05, incident_txt,
    transform=axes[1].transAxes,
    ha='left', va='bottom', fontsize=6.5,
    bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='0.7', lw=0.4))

plt.subplots_adjust(left=0.07, right=0.99, bottom=0.14, top=0.96, wspace=0.12)

out_png = os.path.join(HERE, 'compare_reflionx_xi_scan.png')
out_pdf = os.path.join(HERE, 'compare_reflionx_xi_scan.pdf')
fig.savefig(out_png, dpi=600, bbox_inches='tight')
fig.savefig(out_pdf,            bbox_inches='tight')
with open(os.path.join(HERE, 'plot_inputs.json'), 'w') as f:
    json.dump({'DAO_quantity': 'upper-face outward normal flux 2*pi*sum(w*mu*I), mu>0',
               'normalization': 'Each curve divided by its own mean F_E over 20-50 keV, on the DAO grid.',
               'runs': plot_records}, f, indent=2)
    f.write('\n')
print(f'Saved: {out_png}')
print(f'Saved: {out_pdf}')
if args.show:
    plt.show()
