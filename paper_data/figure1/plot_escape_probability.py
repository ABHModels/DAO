"""Escape probability for subordinate lines: CRD with Voigt wings.

Left:  geometry of the cell-by-cell solver — stacked plasma cells
       bathed by the mean radiation field J_nu, emitting continuum
       j_con directly and line photons after escape factor beta_l.
Right: beta_l(tau_l) = (1 - p_w) beta_K2(tau_l) + p_w
       (Hummer 1982; Ferland et al. 2017).

Author:      Yimin Huang
Affiliation: Fudan University
Email:       huangym23@m.fudan.edu.cn (hyimin0924@gmail.com)
"""

import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

# Resolve output location the same way as paper_data/figure2/plot_kernel_redist.py:
# write next to this script by default, but prefer ./plot/ if it already exists
# (so running from the repo root keeps the historical layout). Works from any cwd.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = "plot" if os.path.isdir("plot") else SCRIPT_DIR
os.makedirs(OUT_DIR, exist_ok=True)
OUT = os.path.join(OUT_DIR, "escape_probability.pdf")
# Palette matched to beamer theme
C_AZURE = "#2A62AF"
C_INK   = "#12203A"
C_VERM  = "#C4302B"
C_GOLD  = "#CDAD6D"
C_SLATE = "#6C7484"
C_MIST  = "#E4E8F0"
C_CELL  = "#CFE0F5"


# ---------------------------------------------------------------
# Cloudy's esca0k2 — Hummer's K_2 one-sided escape probability for
# a pure Doppler profile (a = 0).
#
# Direct line-by-line port of esca0k2() in
#   cloudy/source/rt_escprob.cpp  (lines 425-481)
#
# Rational-approximation coefficients are from
#   Hummer, D.G., 1981, JQSRT, 26, 187
# and the K_2 function itself is the static one-sided escape
# probability that appears as the unmarked lower-envelope curve in
#   Hummer & Rybicki, 1982, ApJ, 254, 767, Fig. 2
# (their Eqs. 2.11-2.12).
#
# Argument convention: input is the LINE-CENTER optical depth taume,
# related to Hummer's tau by tau = taume * sqrt(pi). Limits:
#   taume -> 0   :  esca0k2 -> 1
#   taume >> 1   :  esca0k2 -> 1 / [2 tau sqrt(ln(tau/sqrt(pi)))]
# ---------------------------------------------------------------
SQRTPI = np.sqrt(np.pi)
_A = np.array([1.00, -0.1117897, -0.1249099917, -9.136358767e-3,
               -3.370280896e-4])
_B = np.array([1.00,  0.1566124168, 9.013261660e-3, 1.908481163e-4,
               -1.547417750e-7, -6.657439727e-9])
_C = np.array([1.000, 19.15049608, 100.7986843, 129.5307533, -31.43372468])
_D = np.array([1.00, 19.68910391, 110.2576321, 169.4911399, -16.69969409,
               -36.664480000])


def _esca0k2_scalar(taume):
    tau = taume * SQRTPI
    if tau < 0.01:
        return 1.0 - 2.0 * tau
    elif tau <= 11.0:
        suma = _A[0] + tau*(_A[1] + tau*(_A[2] + tau*(_A[3] + _A[4]*tau)))
        sumb = (_B[0] + tau*(_B[1] + tau*(_B[2] + tau*(_B[3]
                + tau*(_B[4] + _B[5]*tau)))))
        return tau / 2.5066283 * np.log(tau / SQRTPI) + suma / sumb
    else:
        arg = 1.0 / np.log(tau / SQRTPI)
        sumc = _C[0] + arg*(_C[1] + arg*(_C[2] + arg*(_C[3] + _C[4]*arg)))
        sumd = (_D[0] + arg*(_D[1] + arg*(_D[2] + arg*(_D[3]
                + arg*(_D[4] + _D[5]*arg)))))
        return (sumc / sumd) / (2.0 * tau * np.sqrt(np.log(tau / SQRTPI)))


def beta_K2(tau):
    """Hummer K_2 one-sided escape probability (Cloudy esca0k2)."""
    tau = np.atleast_1d(np.asarray(tau, dtype=float))
    return np.array([_esca0k2_scalar(t) for t in tau]).reshape(tau.shape)


# ---------------------------------------------------------------
# Cloudy's esc_PRD_1side — partial (incomplete) redistribution,
# one-sided escape probability with damping wings.
#
# Direct port of esc_PRD_1side() in
#   cloudy/source/rt_escprob.cpp  (lines 116-162)
#
# Form:  beta = 1 / (1 + b(a, tau) * tau),  with
#   atau = a * tau
#   b = 1.6 + 3 (2a)^{-0.12} * f(atau)
#   f(atau) = sqrt(atau)/(1+sqrt(atau))   for atau <= 1
#   f(atau) = 1/(1+atau)                  for atau > 1
#   b capped at 6
# Argument: tau is the line-center optical depth, a is the Voigt
# damping parameter (a = Gamma / (4 pi Delta nu_D)).
# ---------------------------------------------------------------
def _esc_PRD_1side_scalar(tau, a):
    if tau < 0.0:
        return np.nan
    atau = a * tau
    pref = 3.0 * (2.0 * a) ** (-0.12)
    if atau > 1.0:
        b = 1.6 + pref / (1.0 + atau)
    else:
        sa = np.sqrt(atau)
        b = 1.6 + pref * sa / (1.0 + sa)
    b = min(6.0, b)
    return 1.0 / (1.0 + b * tau)


def beta_PRD(tau, a):
    """Cloudy esc_PRD_1side — incomplete-redistribution escape prob."""
    tau = np.atleast_1d(np.asarray(tau, dtype=float))
    return np.array([_esc_PRD_1side_scalar(t, a) for t in tau]).reshape(tau.shape)


fig, (axL, axR) = plt.subplots(
    1, 2, figsize=(10.4, 4.4),
    gridspec_kw={"width_ratios": [1.0, 1.15]},
    constrained_layout=True,
)

# ===============================================================
# Left: celled-plasma geometry
# ===============================================================
axL.set_xlim(0, 1)
axL.set_ylim(0, 1)
axL.set_aspect("equal")
axL.set_axis_off()

# Stack of cells (optical depth grows downward)
cell_tops = [0.88, 0.70, 0.52, 0.34, 0.16]
tau_labels = [r"$\tau_{d-1}$", r"$\tau_{d}$", r"$\tau_{d+1}$",
              r"$\tau_{d+2}$", r"$\tau_{d+3}$"]
x_lo, x_hi = 0.22, 0.86

# tau grid lines + labels
for y, lbl in zip(cell_tops, tau_labels):
    axL.plot([x_lo - 0.03, x_hi + 0.03], [y, y],
             color=C_SLATE, lw=0.8, ls="--", zorder=1)
    axL.text(x_lo - 0.05, y, lbl, ha="right", va="center",
             fontsize=10.5, color=C_INK)

# Fill each cell with plasma (dots)
rng = np.random.default_rng(7)
for i in range(len(cell_tops) - 1):
    yt, yb = cell_tops[i], cell_tops[i + 1]
    highlight = (i == 1)  # cell d is between tau_d and tau_{d+1}
    face = C_CELL if highlight else C_MIST
    axL.add_patch(Rectangle(
        (x_lo, yb), x_hi - x_lo, yt - yb,
        facecolor=face, edgecolor="none", zorder=2, alpha=0.9,
    ))
    # plasma dots
    n_dots = 28
    xs = rng.uniform(x_lo + 0.02, x_hi - 0.02, n_dots)
    ys = rng.uniform(yb + 0.015, yt - 0.015, n_dots)
    axL.scatter(xs, ys, s=7, color=C_SLATE, alpha=0.55, zorder=3)

# Direction-of-tau arrow on the far left
axL.annotate(
    "", xy=(x_lo - 0.15, cell_tops[-1] - 0.02),
    xytext=(x_lo - 0.15, cell_tops[0] + 0.02),
    arrowprops=dict(arrowstyle="->", color=C_SLATE, lw=1.1),
)
axL.text(x_lo - 0.17, (cell_tops[0] + cell_tops[-1]) / 2,
         r"$\tau$", ha="right", va="center",
         fontsize=12, color=C_SLATE)

# ----- Highlighted cell d -----
yt, yb = cell_tops[1], cell_tops[2]   # tau_d to tau_{d+1}
xm = (x_lo + x_hi) / 2
ym = (yt + yb) / 2

# Mean radiation field J̄ — incoming from both sides of the cell
# (gold curved arrows)
for x0 in np.linspace(x_lo + 0.1, x_hi - 0.1, 3):
    axL.add_patch(FancyArrowPatch(
        (x0 - 0.03, yt + 0.06), (x0, yt),
        arrowstyle="-|>", color=C_GOLD, lw=1.3, mutation_scale=10,
        connectionstyle="arc3,rad=0.15", zorder=4,
    ))
    axL.add_patch(FancyArrowPatch(
        (x0 + 0.03, yb - 0.06), (x0, yb),
        arrowstyle="-|>", color=C_GOLD, lw=1.3, mutation_scale=10,
        connectionstyle="arc3,rad=0.15", zorder=4,
    ))

axL.text(x_hi + 0.05, ym + 0.02, r"$\bar J_\nu$",
         color=C_GOLD, fontsize=13, va="center", weight="bold")
axL.text(x_hi + 0.05, ym - 0.04, "(mean\nradiation\nfield)",
         color=C_GOLD, fontsize=7.8, va="top", ha="left")

# Continuum emissivity j_con (blue outgoing arrow, top-left)
axL.add_patch(FancyArrowPatch(
    (x_lo + 0.22, yt), (x_lo + 0.03, yt + 0.12),
    arrowstyle="-|>", color=C_AZURE, lw=2.0, mutation_scale=14,
    connectionstyle="arc3,rad=0.25", zorder=5,
))
axL.text(x_lo - 0.02, yt + 0.13, r"$j_{\rm con}$",
         color=C_AZURE, fontsize=13, weight="bold",
         ha="left", va="bottom")

# Escaped line emissivity β_ℓ j_line (red outgoing arrow, top-right)
axL.add_patch(FancyArrowPatch(
    (x_hi - 0.22, yt), (x_hi - 0.03, yt + 0.12),
    arrowstyle="-|>", color=C_VERM, lw=2.0, mutation_scale=14,
    connectionstyle="arc3,rad=-0.25", zorder=5,
))
axL.text(x_hi - 0.02, yt + 0.13,
         r"$\beta_\ell\, j_{\rm line}$",
         color=C_VERM, fontsize=13, weight="bold",
         ha="right", va="bottom")

# Highlight box outline
axL.add_patch(Rectangle(
    (x_lo, yb), x_hi - x_lo, yt - yb,
    facecolor="none", edgecolor=C_INK, lw=1.6, zorder=6,
))
axL.text(xm, (yt + yb) / 2, "cell $d$",
         ha="center", va="center",
         fontsize=10.5, color=C_INK, style="italic",
         bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                   edgecolor=C_INK, lw=0.7))

# ===============================================================
# Right: beta_l vs tau_l
# ===============================================================
tau = np.logspace(-3, 5, 500)
bK2 = beta_K2(tau)

# Doppler core — solid black
axR.loglog(tau, bK2, color=C_INK, lw=2.0, ls="-",
           label=r"$\beta_{K_2}$")

# CRD + wing-leakage: one curve per p_w, each color explicitly labeled
pw_vals   = [1e-4, 1e-3, 1e-2]
pw_colors = [C_AZURE, C_VERM, C_GOLD]
for pw, c in zip(pw_vals, pw_colors):
    bl = (1.0 - pw) * bK2 + pw
    axR.loglog(tau, bl, color=c, lw=1.8, ls="--",
               label=rf"$\beta_\ell$, $p_w=10^{{{int(np.log10(pw))}}}$")
    axR.axhline(pw, color=c, lw=0.8, ls=":", alpha=0.45)

# Cloudy esc_PRD_1side — single representative curve (a = 1e-3) to
# avoid duplicating colors. Dotted, distinct dark color.
a_prd = 1e-3
bp = beta_PRD(tau, a_prd)
axR.loglog(tau, bp, color="#2A6F4F", lw=1.8, ls=":",
           label=rf"$\beta_{{\rm PRD}}$, $a={a_prd:g}$ ")

axR.set_xlabel(r"Line-center optical depth  $\tau_\ell$", fontsize=11)
axR.set_ylabel(r"Escape probability  $\beta_\ell$", fontsize=11)
axR.set_xlim(1e-3, 1e5)
axR.set_ylim(5e-5, 2.0)
axR.grid(True, which="both", ls=":", alpha=0.4)
axR.tick_params(which="both", direction="in", top=True, right=True)
axR.legend(frameon=False, loc="lower left", fontsize=9.5)

# axR.text(
#     0.7, 0.95,
#     r"$\beta_\ell = (1-p_w)\,\beta_{K_2}(\tau_\ell) + p_w$",
#     transform=axR.transAxes, fontsize=12.5, color=C_INK,
#     ha="center", va="top",
#     bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
#               edgecolor=C_SLATE, lw=0.8),
# )
# axR.text(
#     0.70, 0.83,
#     r"Hummer & Rybicki 1982;",
#     transform=axR.transAxes, fontsize=8.5, color=C_SLATE,
#     ha="center", va="top", style="italic",
# )

plt.show()
fig.savefig(OUT, bbox_inches="tight")
print(f"wrote {OUT}")