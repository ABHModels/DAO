#!/usr/bin/env python
"""Angle-averaged vs directional kernel at fixed incidence (cutoffpl).

Two cutoffpl runs with identical physics (Gamma=2, E_cut=300 keV, nH=10^15,
log xi=2, incidence cos=0.70 -> snapped GL node 0.797). Only the Compton
kernel differs:

  a1f9c7b5  : directional kernel (angsca=true)
  1e89bf78  : angle-averaged kernel (angsca=false), v2

Two-panel figure:
2x2 small multiples, one cell per viewing GL node (mu = 0.183, 0.526,
0.797, 0.960). Each cell has a stacked pair: top panel = raw emergent
I_E for both kernels (solid = directional, dashed = angle-averaged) with
the gap shaded; bottom panel = ratio angle-averaged / directional with a
dotted y=1 line for parity. Band-integrated fractional differences in
1-30 keV and 30-300 keV are annotated on each spectrum panel.
"""

import os
import re
import sys
import glob
import json
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

HASH_DIR = "a1f9c7b5"   # directional reference
HASH_AVG = "1e89bf78"   # angle-averaged (v2)

SPEC_YLIM      = (1e4, 1e14)         # raw I_E range for cutoffpl, log\xi=2
SPEC_XLIM_keV  = (0.05, 500.0)
RATIO_YLIM     = (0.5, 2.0)
BAND_SOFT_keV  = (1.0,  30.0)        # soft band for fractional-diff annotation
BAND_HARD_keV  = (30.0, 300.0)       # hard band for fractional-diff annotation


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


def latest_iter_file(hash_):
    files = glob.glob(f"{HERE}/emergent_{hash_}_iter*.dat")
    if not files:
        sys.exit(f"No emergent_{hash_}_iter*.dat in {HERE}/")
    def idx(f):
        m = re.search(r"_iter(\d+)\.dat$", os.path.basename(f))
        return int(m.group(1)) if m else -1
    f = max(files, key=idx)
    return f, idx(f)


def read_emergent(em_file):
    mus = []
    with open(em_file) as f:
        for line in f:
            if not line.startswith("#"):
                break
            m = re.search(r"I\(mu=([-\d.]+)\)", line)
            if m:
                mus.append(float(m.group(1)))
    mus = np.array(mus)
    data = np.loadtxt(em_file, comments="#")
    E_eV = data[:, 0]
    I = {float(mu): data[:, 3 + i] for i, mu in enumerate(mus)}
    return E_eV, mus, I


def total_flux(E, F):
    m = np.isfinite(F) & (F > 0)
    return float(np.trapz(F[m], E[m])) if m.any() else 1.0


def band_int(E_eV, F, E_lo_keV, E_hi_keV):
    m = (E_eV >= E_lo_keV * 1e3) & (E_eV <= E_hi_keV * 1e3) & np.isfinite(F) & (F > 0)
    return float(np.trapz(F[m], E_eV[m])) if m.any() else 0.0


def main():
    runs = {}
    for label, h in [("directional", HASH_DIR), ("angle-averaged", HASH_AVG)]:
        em, it = latest_iter_file(h)
        E_eV, mus, I = read_emergent(em)
        print(f"  {label:14s}  hash={h}  {os.path.basename(em)} (iter {it})")
        runs[label] = (E_eV, mus, I, it)

    mus_full = runs["directional"][1]
    mu_views = sorted(float(m) for m in mus_full if m > 0)
    print(f"Viewing GL nodes: {[f'{m:.3f}' for m in mu_views]}")

    with open(f"{HERE}/params.json") as f:
        par = json.load(f)

    E_dir, _, I_dir, it_dir = runs["directional"]
    E_avg, _, I_avg, it_avg = runs["angle-averaged"]

    # 2x2 small multiples, each cell is a (spectrum, ratio) stack.
    # Outer GridSpec: 2 rows of cells x 2 cols. Inner per-cell: height_ratios=[3,1].
    fig = plt.figure(figsize=(7.0, 6.4))
    outer = fig.add_gridspec(2, 2, wspace=0.05, hspace=0.18,
                             left=0.09, right=0.985, bottom=0.07, top=0.985)

    col_dir  = "#1f3a93"   # deep blue
    col_avg  = "#d62728"   # red
    col_fill = "#d6c5c5"   # light pink fill between

    sp_axes, hr_axes = [], []
    for k, mu in enumerate(mu_views):
        r, c = divmod(k, 2)
        inner = outer[r, c].subgridspec(2, 1, height_ratios=[3, 1], hspace=0.0)
        ax_sp = fig.add_subplot(inner[0])
        ax_hr = fig.add_subplot(inner[1], sharex=ax_sp)
        sp_axes.append(ax_sp)
        hr_axes.append(ax_hr)

        F_dir = I_dir[mu]
        F_avg = I_avg[mu]
        # interpolate directional onto angle-averaged grid for fill / ratio
        F_dir_on_avg = np.interp(E_avg, E_dir, F_dir)
        E_keV_avg = E_avg / 1.0e3

        # --- spectrum panel ---
        # Clip each curve to the panel floor before filling so the shading
        # still covers the visible gap when one curve dives below SPEC_YLIM[0]
        # at the cutoff (the previous `where` mask dropped the fill there).
        good = np.isfinite(F_dir_on_avg) & np.isfinite(F_avg)
        F_dir_clip = np.maximum(F_dir_on_avg, SPEC_YLIM[0])
        F_avg_clip = np.maximum(F_avg,        SPEC_YLIM[0])
        ax_sp.fill_between(E_keV_avg, F_dir_clip, F_avg_clip, where=good,
                           color=col_fill, alpha=0.7, lw=0, zorder=1)
        ax_sp.loglog(E_dir / 1.0e3, F_dir, color=col_dir, lw=1.1, alpha=0.95,
                     ls='-',  zorder=3, label=f"angle dep")
        ax_sp.loglog(E_keV_avg,    F_avg, color=col_avg, lw=1.1, alpha=0.95,
                     ls='--', zorder=4, label=f"angle avg")
        ax_sp.set_xlim(*SPEC_XLIM_keV)
        ax_sp.set_ylim(*SPEC_YLIM)
        plt.setp(ax_sp.get_xticklabels(), visible=False)

        soft_d = band_int(E_dir, F_dir, *BAND_SOFT_keV)
        soft_a = band_int(E_avg, F_avg, *BAND_SOFT_keV)
        hard_d = band_int(E_dir, F_dir, *BAND_HARD_keV)
        hard_a = band_int(E_avg, F_avg, *BAND_HARD_keV)
        d_soft = (soft_a - soft_d) / soft_d * 100 if soft_d > 0 else float('nan')
        d_hard = (hard_a - hard_d) / hard_d * 100 if hard_d > 0 else float('nan')

        ax_sp.text(0.04, 0.96, fr"viewing $\mu={mu:.3f}$",
                   transform=ax_sp.transAxes, ha="left", va="top",
                   fontsize=7.5, weight='bold')
        # ax_sp.text(0.96, 0.04,
        #            f"1–30 keV:  {d_soft:+.1f}%\n"
        #            f"30–300 keV: {d_hard:+.1f}%",
        #            transform=ax_sp.transAxes, ha="right", va="bottom",
        #            fontsize=6.5, family='monospace',
        #            bbox=dict(boxstyle="round,pad=0.25", fc="white",
        #                      ec="0.7", lw=0.4))

        # --- ratio panel ---
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(good, F_avg / F_dir_on_avg, np.nan)
        ax_hr.plot(E_keV_avg, ratio, color=col_avg, lw=1.0, alpha=0.95)
        ax_hr.axhline(1.0, color="0.4", ls=":", lw=0.6, zorder=0)
        ax_hr.set_ylim(*RATIO_YLIM)
        ax_hr.set_yticks([0.6, 1.0, 1.5])
        ax_hr.set_xscale('log')

        # left-column gets y-labels; bottom-row gets x-label
        if c == 0:
            ax_sp.set_ylabel(r"$I_E$  (erg cm$^{-2}$ s$^{-1}$ eV$^{-1}$ sr$^{-1}$)",
                             fontsize=7.5)
            ax_hr.set_ylabel("avg / dir", fontsize=7)
        else:
            plt.setp(ax_sp.get_yticklabels(), visible=False)
            plt.setp(ax_hr.get_yticklabels(), visible=False)
        if r == 1:
            ax_hr.set_xlabel(r"Energy (keV)")
        else:
            plt.setp(ax_hr.get_xticklabels(), visible=False)

    # single shared legend in the top-left spectrum panel
    sp_axes[0].legend(loc="lower left", borderpad=0.3, fontsize=6.5)

    # parameter box in the top-right spectrum panel (lower-right corner)
    param_txt = (
        f"{par['corona']}\n"
        rf"$\Gamma = {par['Gamma']:g}$" + "\n"
        rf"$E_{{\rm cut}} = {par['E_cut']:g}$ keV" + "\n"
        rf"$n_{{\rm H}} = 10^{{{par['nh']:g}}}$ cm$^{{-3}}$" + "\n"
        rf"$\log\xi = {par['zeta']:g}$" + "\n"
        rf"$\mu_{{\rm inc}} = 0.797$"
    )
    sp_axes[1].text(0.04, 0.04, param_txt, transform=sp_axes[1].transAxes,
                    ha="left", va="bottom", fontsize=6.5,
                    bbox=dict(boxstyle="round,pad=0.3", fc="white",
                              ec="0.7", lw=0.4))

    out_png = f"{HERE}/angavg_vs_directional.png"
    out_pdf = f"{HERE}/angavg_vs_directional.pdf"
    fig.savefig(out_png, dpi=600, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    print(f"Saved: {out_png}")
    print(f"Saved: {out_pdf}")


if __name__ == "__main__":
    main()
