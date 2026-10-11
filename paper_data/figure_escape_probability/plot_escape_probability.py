"""Revised escape-probability figure.

This script reconstructs the later paper-style version of the line escape
figure.  The important semantic change from the original plotting script is
that continuum destruction and collisional de-excitation are not plotted as
positive contributors to the escaping/surviving line power.  They enter only
through the denominator of the effective survival probability,

    P_ul = (beta + P_elec) (1 + y) / (beta + P_elec + P_dest + y).

The bottom panels therefore decompose P_ul only into the two surviving
channels:

    P_beta = beta   (1 + y) / D,
    P_elec = P_elec (1 + y) / D,

where D = beta + P_elec + P_dest + y.

The slab-averaged value printed in each top panel is width-weighted in Thomson
depth, because the tabulated line powers are local power densities.  The
selected diagnostic file stores cell-midpoint tau_T values, so the cell edges
are reconstructed from tau_mid = (tau_left + tau_right) / 2 with tau_left = 0
at the illuminated surface.
"""

from __future__ import annotations

import os
import re
import shutil
from dataclasses import dataclass

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SELECTED_FILE = os.path.join(HERE, "line_escape_selected_64e9e118_iter011.dat")
OUT_PDF = os.path.join(HERE, "escape_probability.pdf")
OUT_PNG = os.path.join(HERE, "escape_probability.png")
DRAFT_FIGURE_DIR = os.path.abspath(os.path.join(HERE, "..", "..", "Figure"))

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import LogLocator, NullFormatter


@dataclass
class LineProfile:
    selected: int
    ip: int
    energy_ev: float
    raw_label: str
    comment: str
    tau: np.ndarray
    thin: np.ndarray
    escaped: np.ndarray
    beta: np.ndarray
    pelec: np.ndarray
    pdest: np.ndarray
    y: np.ndarray


ROMAN = [
    "",
    "I",
    "II",
    "III",
    "IV",
    "V",
    "VI",
    "VII",
    "VIII",
    "IX",
    "X",
    "XI",
    "XII",
    "XIII",
    "XIV",
    "XV",
    "XVI",
    "XVII",
    "XVIII",
    "XIX",
    "XX",
    "XXI",
    "XXII",
    "XXIII",
    "XXIV",
    "XXV",
    "XXVI",
    "XXVII",
]

C_THIN = "#56585c"
C_SURVIVE = "#2a78d6"
C_REMOVED = "#969ca6"
C_TOTAL = "#111417"
C_LINE = "#2a78d6"
C_ELEC = "#2a9d75"
C_GRID = "#d8dce2"
TINY = 1.0e-6


def read_selected_profiles(path: str) -> list[LineProfile]:
    groups: dict[int, list[dict[str, object]]] = {}
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 14:
                continue
            selected = int(parts[0])
            groups.setdefault(selected, []).append(
                {
                    "depth": int(parts[1]),
                    "tau": float(parts[2]),
                    "ip": int(parts[3]),
                    "energy_ev": float(parts[4]),
                    "thin": float(parts[5]),
                    "escaped": float(parts[6]),
                    "beta": float(parts[9]),
                    "pelec": float(parts[10]),
                    "pdest": float(parts[11]),
                    "y": float(parts[12]),
                    "raw_label": parts[13],
                    "comment": parts[14] if len(parts) > 14 else "",
                }
            )

    profiles: list[LineProfile] = []
    for selected in sorted(groups):
        rows = sorted(groups[selected], key=lambda row: int(row["depth"]))
        profiles.append(
            LineProfile(
                selected=selected,
                ip=int(rows[0]["ip"]),
                energy_ev=float(rows[0]["energy_ev"]),
                raw_label=str(rows[0]["raw_label"]),
                comment=str(rows[0]["comment"]),
                tau=np.array([float(row["tau"]) for row in rows]),
                thin=np.array([float(row["thin"]) for row in rows]),
                escaped=np.array([float(row["escaped"]) for row in rows]),
                beta=np.array([float(row["beta"]) for row in rows]),
                pelec=np.array([float(row["pelec"]) for row in rows]),
                pdest=np.array([float(row["pdest"]) for row in rows]),
                y=np.array([float(row["y"]) for row in rows]),
            )
        )
    if len(profiles) < 4:
        raise RuntimeError(f"expected at least four selected profiles in {path}")
    return profiles[:4]


def tau_cell_widths_from_midpoints(tau_mid: np.ndarray) -> np.ndarray:
    """Recover d tau from midpoint coordinates.

    The diagnostic stores tau_mid rather than explicit cell edges.  Assuming
    tau_mid[i] = (edge[i] + edge[i + 1]) / 2 and edge[0] = 0 reconstructs the
    grid used for the Sep. 2026 figure averages.
    """

    if tau_mid.ndim != 1 or tau_mid.size == 0:
        raise ValueError("tau_mid must be a non-empty 1D array")
    edges = np.zeros(tau_mid.size + 1)
    for i, midpoint in enumerate(tau_mid):
        edges[i + 1] = 2.0 * midpoint - edges[i]
    widths = np.diff(edges)
    if np.any(widths <= 0.0):
        raise ValueError("reconstructed non-positive tau cell width")
    return widths


def spectroscopic_label(raw_label: str, comment: str) -> str:
    match = re.match(r"\s*([A-Z][a-z]?)\s*(\d+)\s+([\d.]+)\s*([A-Za-z]*)", raw_label)
    if not match:
        return " ".join(raw_label.split())

    elem = match.group(1)
    stage = int(match.group(2))
    wavelength = float(match.group(3))
    unit = match.group(4)
    roman = ROMAN[stage] if 0 < stage < len(ROMAN) else str(stage)
    unit_text = r"$\AA$" if unit.upper().startswith("A") else unit

    tag = ""
    if "H-like" in comment and "1^2S" in comment:
        if "2^2P" in comment:
            tag = r" Ly$\alpha$"
        elif "n=  3" in comment:
            tag = r" Ly$\beta$"

    if wavelength < 10.0:
        wavelength_text = f"{wavelength:.4f}"
    else:
        wavelength_text = f"{wavelength:.2f}" if wavelength < 100.0 else f"{wavelength:.1f}"
    return f"{elem} {roman}{tag} {wavelength_text} {unit_text}"


def survival_components(profile: LineProfile) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    den = profile.beta + profile.pelec + profile.pdest + profile.y
    total = np.zeros_like(den)
    line_part = np.zeros_like(den)
    elec_part = np.zeros_like(den)
    valid = den > 0.0
    boost = 1.0 + profile.y[valid]
    total[valid] = (profile.beta[valid] + profile.pelec[valid]) * boost / den[valid]
    line_part[valid] = profile.beta[valid] * boost / den[valid]
    elec_part[valid] = profile.pelec[valid] * boost / den[valid]
    return total, line_part, elec_part


def width_weighted_pbar(profile: LineProfile) -> float:
    widths = tau_cell_widths_from_midpoints(profile.tau)
    denominator = np.sum(profile.thin * widths)
    if denominator <= 0.0:
        return np.nan
    return float(np.sum(profile.escaped * widths) / denominator)


def pbar_text(value: float) -> str:
    if not np.isfinite(value):
        return r"$\overline{P}_{ul}=\cdots$"
    if value < 0.01:
        mantissa, exponent = f"{value:.2e}".split("e")
        return rf"$\overline{{P}}_{{ul}}={float(mantissa):.2f}\times10^{{{int(exponent)}}}$"
    return rf"$\overline{{P}}_{{ul}}={value:.3f}$"


def main() -> None:
    profiles = read_selected_profiles(SELECTED_FILE)

    plt.rcParams.update(
        {
            "font.size": 7.4,
            "axes.linewidth": 0.7,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "mathtext.fontset": "cm",
            "axes.labelpad": 2.0,
            "legend.handlelength": 1.7,
        }
    )

    fig, axes = plt.subplots(
        2,
        4,
        figsize=(7.23, 4.756),
        sharex=True,
        sharey="row",
        gridspec_kw={"height_ratios": [1.0, 1.0], "hspace": 0.11, "wspace": 0.08},
    )
    fig.subplots_adjust(left=0.082, right=0.988, top=0.875, bottom=0.155)

    letters = ["(a)", "(b)", "(c)", "(d)"]
    for col, profile in enumerate(profiles):
        tau = profile.tau
        thin = profile.thin
        escaped = profile.escaped
        emitting = thin > 0.0
        surviving = escaped > 0.0

        ax_top = axes[0, col]
        both = emitting & surviving
        ax_top.fill_between(
            tau[both],
            np.maximum(escaped[both], TINY),
            thin[both],
            where=thin[both] >= escaped[both],
            color=C_REMOVED,
            alpha=0.24,
            linewidth=0.0,
            zorder=1,
        )
        ax_top.plot(
            tau[emitting],
            thin[emitting],
            color=C_THIN,
            linewidth=1.15,
            linestyle=(0, (4, 2)),
            zorder=3,
        )
        ax_top.plot(
            tau[surviving],
            escaped[surviving],
            color=C_SURVIVE,
            linewidth=1.7,
            zorder=4,
        )
        ax_top.set_xscale("log")
        ax_top.set_yscale("log")
        ax_top.set_title(
            f"{letters[col]} {spectroscopic_label(profile.raw_label, profile.comment)}",
            loc="left",
            fontsize=7.6,
            pad=4.0,
        )
        ax_top.text(
            0.05,
            0.065,
            pbar_text(width_weighted_pbar(profile)),
            transform=ax_top.transAxes,
            fontsize=7.4,
            color=C_TOTAL,
            ha="left",
            va="bottom",
        )

        total, line_part, elec_part = survival_components(profile)
        ax_bot = axes[1, col]
        mask = emitting
        ax_bot.plot(
            tau[mask],
            np.clip(total[mask], TINY, 2.0),
            color=C_TOTAL,
            linewidth=1.75,
            zorder=4,
        )
        ax_bot.plot(
            tau[mask],
            np.clip(line_part[mask], TINY, 2.0),
            color=C_LINE,
            linewidth=1.2,
            linestyle=(0, (4, 2)),
            zorder=6,
        )
        ax_bot.plot(
            tau[mask],
            np.clip(elec_part[mask], TINY, 2.0),
            color=C_ELEC,
            linewidth=1.2,
            linestyle=(0, (4, 2, 1.2, 2)),
            zorder=5,
        )
        ax_bot.axhline(1.0, color=C_GRID, linewidth=0.65, linestyle=(0, (1, 2)), zorder=1)
        ax_bot.set_xscale("log")
        ax_bot.set_yscale("log")
        ax_bot.set_ylim(TINY, 2.0)

        for ax in (ax_top, ax_bot):
            ax.grid(True, which="major", color=C_GRID, linewidth=0.5, zorder=0)
            ax.tick_params(which="both", top=True, right=True, labelsize=6.6, pad=2.0)
            ax.xaxis.set_major_locator(LogLocator(numticks=5))
            ax.xaxis.set_minor_locator(LogLocator(subs=np.arange(2, 10) * 0.1, numticks=100))
            ax.xaxis.set_minor_formatter(NullFormatter())

    tau_positive = np.concatenate([p.tau[p.tau > 0.0] for p in profiles])
    axes[0, 0].set_xlim(tau_positive.min() * 0.75, tau_positive.max() * 1.12)

    all_thin = np.concatenate([p.thin[p.thin > 0.0] for p in profiles])
    all_escaped = np.concatenate([p.escaped[p.escaped > 0.0] for p in profiles])
    axes[0, 0].set_ylim(
        max(all_escaped.min() * 0.45, all_thin.max() * 1.0e-6),
        all_thin.max() * 3.0,
    )

    axes[0, 0].set_ylabel(
        r"Local line-power density" + "\n" + r"[erg cm$^{-3}$ s$^{-1}$]",
        fontsize=7.6,
    )
    axes[1, 0].set_ylabel(
        r"Fraction of optically" + "\n" + r"thin line power",
        fontsize=7.6,
    )
    for ax in axes[1, :]:
        ax.set_xlabel(r"Thomson depth $\tau_{\rm T}$", fontsize=7.4)

    top_handles = [
        Line2D([0], [0], color=C_THIN, linewidth=1.15, linestyle=(0, (4, 2)),
               label="optically thin line power"),
        Line2D([0], [0], color=C_SURVIVE, linewidth=1.7,
               label="passed to continuum/Compton solver"),
        Patch(facecolor=C_REMOVED, alpha=0.24,
              label="removed before continuum transfer"),
    ]
    fig.legend(
        handles=top_handles,
        ncol=3,
        frameon=False,
        fontsize=7.0,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.982),
        columnspacing=1.05,
        handletextpad=0.45,
    )

    bottom_handles = [
        Line2D([0], [0], color=C_TOTAL, linewidth=1.75,
               label=r"$P_{ul}$ total surviving fraction"),
        Line2D([0], [0], color=C_LINE, linewidth=1.2, linestyle=(0, (4, 2)),
               label="line-escape part"),
        Line2D([0], [0], color=C_ELEC, linewidth=1.2, linestyle=(0, (4, 2, 1.2, 2)),
               label="electron-scattering part"),
    ]
    fig.legend(
        handles=bottom_handles,
        ncol=3,
        frameon=False,
        fontsize=7.0,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.026),
        columnspacing=1.25,
        handletextpad=0.45,
    )

    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PNG, dpi=240)
    os.makedirs(DRAFT_FIGURE_DIR, exist_ok=True)
    shutil.copyfile(OUT_PDF, os.path.join(DRAFT_FIGURE_DIR, os.path.basename(OUT_PDF)))
    shutil.copyfile(OUT_PNG, os.path.join(DRAFT_FIGURE_DIR, os.path.basename(OUT_PNG)))

    print(f"wrote {OUT_PDF}")
    print(f"wrote {OUT_PNG}")
    print(f"synced {DRAFT_FIGURE_DIR}")
    print(f"selected {SELECTED_FILE}")
    for profile in profiles:
        print(
            f"selected={profile.selected} ip={profile.ip} "
            f"{spectroscopic_label(profile.raw_label, profile.comment)} "
            f"Pbar_width_weighted={width_weighted_pbar(profile):.6e}"
        )


if __name__ == "__main__":
    main()
