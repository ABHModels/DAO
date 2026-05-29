"""
Compton redistribution kernel K(theta_out, E_out) at fixed
(E_in, mu_in, T), evaluated via compton_kernel_element().

Two panels, shared scales (axes + color), highlight the energy
bands the reader cares about:

  (a) Iron K-alpha line:        E_in = 6.4   keV
  (b) Compton-hump source:      E_in = 100   keV

Both panels share the same E_out range (1-500 keV, log scale)
and the same color scale, so the relative magnitudes can be
read off directly.

Reads the two .dat slices that sit next to this script and writes the
figure into the same folder. No arguments; works from any cwd:
    python plot_kernel_redist.py

Author:      Yimin Huang
Affiliation: Fudan University
Email:       huangym23@m.fudan.edu.cn (hyimin0924@gmail.com)
"""

import os
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import LogNorm

# ------------------------------------------------------------------
# Nature-style typography
# ------------------------------------------------------------------
plt.rcParams.update({
    "font.family"      : "sans-serif",
    "font.sans-serif"  : ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size"        : 8,
    "axes.titlesize"   : 9,
    "axes.labelsize"   : 8,
    "axes.linewidth"   : 0.8,
    "xtick.labelsize"  : 7.5,
    "ytick.labelsize"  : 7.5,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "xtick.major.size" : 3.0,
    "ytick.major.size" : 3.0,
    "legend.fontsize"  : 7,
    "legend.frameon"   : False,
    "savefig.dpi"      : 600,
    "figure.dpi"       : 200,
    "mathtext.fontset" : "stixsans",
})

# Compton-hump energy band (eV)
HUMP_BAND  = (20e3, 40e3)
FE_BAND    = (6.20e3, 6.45e3)


def read_slice(path):
    meta, data_rows = {}, []
    mu_out = Eout = None
    with open(path) as f:
        for ln in f:
            if ln.startswith("#"):
                m = re.match(r"#\s*(\w+)\s*=\s*(.*)", ln)
                if not m: continue
                key, val = m.group(1), m.group(2).strip()
                if   key == "mu_out_grid":  mu_out = np.array([float(x) for x in val.split()])
                elif key == "Eout_grid_eV": Eout   = np.array([float(x) for x in val.split()])
                else:                       meta[key] = val
            elif ln.strip():
                data_rows.append([float(x) for x in ln.split()])
    K = np.array(data_rows)
    return dict(
        T     = float(meta["T"].split()[0]),
        E_in  = float(meta["E_in"].split()[0]),
        mu_in = float(meta["mu_in"]),
        mu_out= mu_out,
        Eout  = Eout,
        K     = K,
    )


def edges(c):
    e = np.empty(len(c) + 1)
    e[1:-1] = 0.5 * (c[1:] + c[:-1])
    e[0]    = c[0]  - 0.5 * (c[1]  - c[0])
    e[-1]   = c[-1] + 0.5 * (c[-1] - c[-2])
    return e


def draw_panel(ax, slc, vmin, vmax, band, band_label, mark_upscatter=False):
    theta_deg = np.degrees(np.arccos(slc["mu_out"]))
    order = np.argsort(theta_deg)
    theta_sorted = theta_deg[order]
    K_sorted = slc["K"][order, :]

    theta_edges = edges(theta_sorted)
    Eout_edges  = 10 ** edges(np.log10(slc["Eout"]))

    K_plot = np.where(K_sorted > 0, K_sorted, np.nan)

    cmap = plt.cm.magma.copy()
    cmap.set_bad(cmap(0.0))      # NaN -> bottom of colormap (no white blocks)

    pcm = ax.pcolormesh(
        theta_edges, Eout_edges, K_plot.T,
        norm=LogNorm(vmin=vmin, vmax=vmax),
        cmap=cmap, shading="auto", rasterized=True,
    )

    # Energy band
    ax.axhspan(band[0], band[1], color="#ffd166", alpha=0.22, lw=0)
    ax.axhline(band[0], color="#ffd166", lw=0.7, alpha=0.85)
    ax.axhline(band[1], color="#ffd166", lw=0.7, alpha=0.85)
    ax.text(176, np.sqrt(band[0]*band[1]), band_label,
            color="#ffe08a", fontsize=9.5, fontweight="bold",
            ha="right", va="center",
            path_effects=[pe.withStroke(linewidth=2.0, foreground="black")])

    # E_in marker
    ax.axhline(slc["E_in"], color="white", ls=":", lw=0.8, alpha=0.85)
    ax.text(2, slc["E_in"] * 1.08,
            fr"$E_{{\rm in}}={slc['E_in']/1e3:g}$ keV",
            color="white", fontsize=7.5, ha="left", va="bottom")

    # Incoming-direction marker
    th_in = np.degrees(np.arccos(slc["mu_in"]))
    ax.axvline(th_in, color="white", ls=":", lw=0.8, alpha=0.7)

    # Optional: highlight detailed-balance up-scatter excess at large theta
    if mark_upscatter:
        x0 = 130
        ax.annotate("",
                    xy=(x0, slc["E_in"] * 2.2), xytext=(x0, slc["E_in"]),
                    arrowprops=dict(arrowstyle="-|>", color="#7ee0ff",
                                    lw=1.4, mutation_scale=12))
        ax.annotate("",
                    xy=(x0, slc["E_in"] / 2.2), xytext=(x0, slc["E_in"]),
                    arrowprops=dict(arrowstyle="-|>", color="#7ee0ff",
                                    lw=1.0, mutation_scale=10, alpha=0.55))
        ax.text(x0 + 4, slc["E_in"] * 2.6, "up-scatter\n(favored)",
                color="#7ee0ff", fontsize=8.5, fontweight="bold",
                ha="left", va="bottom",
                path_effects=[pe.withStroke(linewidth=2.0, foreground="black")])
        ax.text(x0 + 4, slc["E_in"] / 2.6, "down-scatter",
                color="#7ee0ff", fontsize=7.5, alpha=0.85,
                ha="left", va="top",
                path_effects=[pe.withStroke(linewidth=1.6, foreground="black")])

    ax.set_yscale("log")
    ax.set_xlim(0, 180)
    ax.set_xticks([0, 30, 60, 90, 120, 150, 180])
    ax.set_xticklabels(["0°", "30°", "60°", "90°",
                        "120°", "150°", "180°"])
    ax.set_xlabel(r"outgoing direction  $\theta_{\rm out}=\arccos\mu_{\rm out}$")
    return pcm


# ------------------------------------------------------------------
# Inputs live next to this script, so paths resolve from any cwd.
# ------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

fe   = read_slice(os.path.join(SCRIPT_DIR, "kernel_redist2d_FeKa.dat"))
hump = read_slice(os.path.join(SCRIPT_DIR, "kernel_redist2d_hump.dat"))

K_all = np.concatenate([fe["K"][fe["K"] > 0], hump["K"][hump["K"] > 0]])
vmax  = K_all.max()
vmin  = vmax * 1e-6

# ------------------------------------------------------------------
# Figure (Nature double-column ~ 183 mm = 7.2 in)
# ------------------------------------------------------------------
fig, axes = plt.subplots(
    1, 2, figsize=(7.2, 3.4), sharey=True,
    gridspec_kw={"wspace": 0.05},
)

pcm = draw_panel(axes[0], fe,   vmin, vmax,
                 band=FE_BAND,   band_label=r"Fe K$\alpha$")
draw_panel(axes[1], hump, vmin, vmax,
           band=HUMP_BAND, band_label="Compton hump",
           mark_upscatter=True)

axes[0].set_ylabel(r"scattered energy $E_{\rm out}$  [eV]")
axes[0].set_ylim(1e3, 5e5)

# Bold panel labels (a), (b)  — Nature style
for ax, lbl in zip(axes, ["a", "b"]):
    ax.text(-0.10, 1.02, lbl, transform=ax.transAxes,
            fontsize=11, fontweight="bold", va="bottom", ha="left")

axes[0].set_title(r"$E_{\rm in}=6.4$ keV   (Fe K$\alpha$)", pad=4)
axes[1].set_title(r"$E_{\rm in}=40$ keV   (Compton-hump band)", pad=4)

# Single shared colorbar on the right
cax = fig.add_axes([0.92, 0.18, 0.013, 0.70])
cbar = fig.colorbar(pcm, cax=cax)
cbar.set_label(r"$K\;[\,1/(m_e c^2)\,]$", fontsize=8)
cbar.ax.tick_params(labelsize=7)

# Layout
fig.subplots_adjust(left=0.07, right=0.91, bottom=0.16, top=0.92)

# Write the figure next to this script.
out = os.path.join(SCRIPT_DIR, "kernel_redist_FeKa_hump.png")
plt.show()
fig.savefig(out, bbox_inches="tight")
fig.savefig(out.replace(".png", ".pdf"), bbox_inches="tight")
print(f"Saved {out} (+ pdf)")
