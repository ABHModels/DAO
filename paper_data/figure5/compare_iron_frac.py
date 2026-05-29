"""Figure 5 — Iron ion-fraction colormap vs Thomson optical depth.

Two panels, iron ionization fraction (Fe I … Fe XXVII) vs Thomson optical
depth τ_T for two versions of the DAO model:

  (a) DAOv1.0 — Fe fractions from abund_9bd2d8c227_C3.fits, τ from the last
      200 rows of temp_9bd2d8c227_C3.dat.
  (b) DAOv2.0 — Fe fractions from the latest iteration of the per-depth
      c63d0c48<N>.iron files, τ from profile_c63d0c48_iter018.dat.

Self-contained: every input file sits next to this script and the figure is
written into the same folder. No arguments:
    python compare_iron_frac.py

Author:      Yimin Huang
Affiliation: Fudan University
Email:       huangym23@m.fudan.edu.cn (hyimin0924@gmail.com)
"""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.ticker import FixedLocator, FixedFormatter
from matplotlib.colors import LogNorm
from matplotlib.cm import ScalarMappable

# ── Nature-style rcParams (matches plot_compare_reflionx.py) ─────────────────
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
    'pdf.fonttype':       42,
    'ps.fonttype':        42,
})

# Self-contained: all input files sit next to this script.
HERE = os.path.dirname(os.path.abspath(__file__))
DAO2_HASH = 'c63d0c48'   # DAOv2.0 run (per-depth .iron files + profile)

# Colormap with zero/masked cells (frac=0, masked by LogNorm) painted with the
# colormap's bottom colour, so they read as "negligible" instead of blank white.
CMAP = plt.get_cmap('cividis').copy()
CMAP.set_bad(CMAP(0.0))
CMAP.set_under(CMAP(0.0))
NORM = LogNorm(vmin=1e-4, vmax=1.0)


def to_roman(num):
    val_syms = [(10, "X"), (9, "IX"), (5, "V"), (4, "IV"), (1, "I")]
    roman = "X" * (num // 10)
    num %= 10
    for val, sym in val_syms:
        while num >= val:
            roman += sym
            num -= val
    return roman


def tau_boundaries(tau):
    """Cell boundaries (in log space) for pcolormesh flat shading."""
    log_tau = np.log10(tau)
    mid = (log_tau[:-1] + log_tau[1:]) / 2.0
    first = log_tau[0] - (mid[0] - log_tau[0])
    last = log_tau[-1] + (log_tau[-1] - mid[-1])
    return 10 ** np.concatenate([[first], mid, [last]])


# ── Model A: DAOv1.0 ─────────────────────────────────────────────────────────
def load_C3():
    tau = pd.read_csv(os.path.join(HERE, 'temp_9bd2d8c227_C3.dat'),
                      sep=r"\s+", header=None).iloc[-200:][0].to_numpy(dtype=float)
    abund = fits.open(os.path.join(HERE, 'abund_9bd2d8c227_C3.fits'))[1].data
    names = [n for n in abund.columns.names if n.startswith('fe_')]
    frac = np.array([abund[n] for n in names], dtype=float)   # [nion, ndepth]
    frac /= frac.sum(axis=0)
    return tau, frac


# ── Model B: DAOv2.0 (one .iron file per depth, last block = latest iter) ────
def load_c63():
    prof = np.loadtxt(os.path.join(HERE, f'profile_{DAO2_HASH}_iter018.dat'))
    tau = prof[:, 1]                                          # tau_mid [Thomson]
    n_depth = prof.shape[0]

    rows = []
    for i in range(n_depth):
        path = os.path.join(HERE, f'{DAO2_HASH}{i}.iron')
        data_lines = [ln for ln in open(path)
                      if ln.strip() and not ln.lstrip().startswith('#')]
        vals = np.array(data_lines[-1].split(), dtype=float)  # latest iteration
        rows.append(vals[1:])                                 # drop depth column
    frac = np.array(rows).T                                   # [nion, ndepth]
    frac /= frac.sum(axis=0)
    return tau, frac


def plot_panel(ax, tau, frac, title, panel_label):
    ion_levels = np.arange(1, frac.shape[0] + 2)
    im = ax.pcolormesh(tau_boundaries(tau), ion_levels, frac,
                       shading='flat', cmap=CMAP, norm=NORM)
    ax.set_xscale('log')
    ax.set_xlabel(r'$\tau_{T}$')
    ax.set_title(title)

    major_ions = [1, 5, 10, 15, 20, 26]
    ax.yaxis.set_major_locator(FixedLocator([i + 0.5 for i in major_ions]))
    ax.yaxis.set_major_formatter(
        FixedFormatter([f"Fe {to_roman(i)}" for i in major_ions]))
    ax.yaxis.set_minor_locator(
        FixedLocator([i + 0.5 for i in range(1, frac.shape[0] + 1)]))
    ax.tick_params(axis='both', which='both')

    ax.text(0.03, 0.96, panel_label, transform=ax.transAxes,
            fontsize=9, fontweight='bold', va='top', ha='left', color='w')
    return im


tau_a, frac_a = load_C3()
tau_b, frac_b = load_c63()

fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), layout="constrained")
fig.supylabel(r'Ionization state')

im = plot_panel(axes[0], tau_a, frac_a, 'DAOv1.0', 'a')
plot_panel(axes[1], tau_b, frac_b, 'DAOv2.0', 'b')

sm = ScalarMappable(cmap=CMAP, norm=NORM)
sm.set_array([])
cbar = fig.colorbar(sm, ax=axes, orientation='horizontal', aspect=40, pad=0.02)
cbar.set_label('Ionization fraction')

plt.show()
out_png = os.path.join(HERE, 'compare_iron_frac.png')
out_pdf = os.path.join(HERE, 'compare_iron_frac.pdf')
fig.savefig(out_png, dpi=600, bbox_inches='tight')
fig.savefig(out_pdf, bbox_inches='tight')
print('Saved:', out_png)
print('Saved:', out_pdf)
