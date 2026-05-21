"""Figure 3 — DAO vs reflionx vs xillvercp over a log xi scan.

Author:      Yimin Huang
Affiliation: Fudan University
Email:       huangym23@m.fudan.edu.cn (hyimin0924@gmail.com)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import re
import json
import glob
import os

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

# ── Gauss–Legendre weights for emergent flux integration ─────────────────────
def gl_weights(nodes):
    n = len(nodes)
    w = np.zeros(n)
    for i, x in enumerate(nodes):
        p0, p1 = 1.0, x
        for j in range(2, n + 1):
            p2 = ((2*j - 1) * x * p1 - (j - 1) * p0) / j
            p0, p1 = p1, p2
        dp = n * (x * p1 - p0) / (x * x - 1.0)
        w[i] = 2.0 / ((1.0 - x * x) * dp * dp)
    return w


def load_emergent_flux(em_file):
    """Return (E [eV], F_up [erg/cm^2/s/eV]) integrated over outgoing mu>0."""
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
    mu_a = np.array(mu_v)
    wt_a = gl_weights(mu_a)
    F = np.zeros_like(E)
    for i in range(len(mu_a)):
        if mu_a[i] > 0:
            F += wt_a[i] * data[:, 3 + i]
    F *= 0.5  
    return E, F


def latest_iter_file(hash_):
    """Return the emergent_iter*.dat file with the highest iteration index."""
    pattern = f'/Users/ym.huang/cloudy_test/results/{hash_}/emergent_iter*.dat'
    files = glob.glob(pattern)
    if not files:
        raise FileNotFoundError(f'No emergent_iter*.dat in {hash_}')
    def idx(f):
        m = re.search(r'emergent_iter(\d+)\.dat$', os.path.basename(f))
        return int(m.group(1)) if m else -1
    return max(files, key=idx)


def load_run(hash_):
    """Load DAO emergent spectrum (latest iteration), reflionx, xillver, params."""
    em_file = latest_iter_file(hash_)
    print(f'[{hash_}] using {os.path.basename(em_file)}')
    E_eV, F_up = load_emergent_flux(em_file)

    with open(f'/Users/ym.huang/cloudy_test/results/{hash_}/params.json') as f:
        par = json.load(f)

    xv = np.loadtxt(f'/Users/ym.huang/cloudy_test/benchmark/xillver/{hash_}.dat')
    E_xv = xv[:, 0] * 1e3                    # keV -> eV
    EFE_xv = xv[:, 2] / E_xv                 # column 2 is E·F_E → divide for F_E

    rf = np.loadtxt(f'/Users/ym.huang/cloudy_test/benchmark/reflionx/'
                    f'spectra_{hash_}.dat')
    E_rf = rf[:, 0] * 1e3
    EFE_rf = rf[:, 1] / E_rf

    return dict(par=par, E=E_eV, F=F_up,
                E_xv=E_xv, F_xv=EFE_xv,
                E_rf=E_rf, F_rf=EFE_rf)


def interp_log(E_target, E_src, F_src):
    return np.exp(np.interp(np.log(E_target),
                            np.log(E_src),
                            np.log(np.clip(F_src, 1e-300, None))))


# ── Runs to plot, ordered by increasing log ξ ────────────────────────────────
HASHES = ['c63d0c48', '8b877447', '50bd6dd6']
PANEL_LABELS = ['a', 'b', 'c']
runs = [load_run(h) for h in HASHES]

# Sort by log xi to guarantee monotonic ordering
order = np.argsort([r['par']['zeta'] for r in runs])
runs = [runs[i] for i in order]
HASHES = [HASHES[i] for i in order]

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
    mean_dao = np.trapz(F[band],    E[band]) / dE_band
    mean_xv  = np.trapz(F_xv[band], E[band]) / dE_band
    mean_rf  = np.trapz(F_rf[band], E[band]) / dE_band

    F_n     = F    / mean_dao
    F_xv_n  = F_xv / mean_xv
    F_rf_n  = F_rf / mean_rf

    # Shade the normalisation band
    ax.axvspan(E_LO, E_HI, color='0.85', alpha=0.5, lw=0, zorder=0)

    ax.loglog(E, F_n,    color=C_DAO,  lw=0.55, alpha=0.85,
              label='DAOv2.0', zorder=3)
    ax.loglog(E, F_rf_n, color=C_REFL, lw=0.9,  alpha=0.95,
              label='reflionx', zorder=5)
    ax.loglog(E, F_xv_n, color=C_XILL, lw=0.9,  alpha=0.95, ls='--',
              label='xillvercp',  zorder=5)

    ax.set_xlim(1e2, 1e6)
    ax.set_ylim(1e-5, 1e3)
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

# In-figure normalisation note (bottom-right of last panel, always visible)
axes[-1].text(
    0.97, 0.05,
    r'$\hat{F}_E = F_E \, / \, \langle F_E\rangle_{20-50\,\mathrm{keV}}$',
    transform=axes[-1].transAxes,
    ha='right', va='bottom', fontsize=6.5,
    bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='0.7', lw=0.4))

# Incident-spectrum parameters as in-figure text (panel a, bottom-left)
par0 = runs[0]['par']
incident_txt = (
    'nthcomp\n'
    rf'$\Gamma = {par0["Gamma"]}$' + '\n'
    rf'$kT_e = {par0["kT_e"]}\,$keV' + '\n'
    rf'$kT_{{\rm bb}} = {par0["kT_bb"]}\,$keV' + '\n'
    rf'$n_{{\rm H}} = 10^{{{par0["nh"]}}}\,$cm$^{{-3}}$' + '\n'
    rf'$\mu_{{\rm inc}} = {{\rm cos45}}^{{\circ}}$'
)
axes[1].text(
    0.03, 0.05, incident_txt,
    transform=axes[1].transAxes,
    ha='left', va='bottom', fontsize=6.5,
    bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='0.7', lw=0.4))

plt.subplots_adjust(left=0.07, right=0.99, bottom=0.14, top=0.96, wspace=0.06)

plt.show()
out_png = '/Users/ym.huang/cloudy_test/image/compare_reflionx_xi_scan.png'
out_pdf = '/Users/ym.huang/cloudy_test/image/compare_reflionx_xi_scan.pdf'
fig.savefig(out_png, dpi=600, bbox_inches='tight')
fig.savefig(out_pdf,            bbox_inches='tight')
print(f'Saved: {out_png}')
print(f'Saved: {out_pdf}')
