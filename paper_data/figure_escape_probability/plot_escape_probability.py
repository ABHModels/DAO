"""Line escape in DAO: how the escape-probability branching reshapes line power.

Two-row small-multiple figure (4 columns = 4 representative lines chosen by
write_line_escape_diagnostics(), one per escape-physics regime):

  Top row    -- WHAT it does: local line power vs Thomson depth, comparing the
                 optically-thin emissivity with the power that actually escapes;
                 the shaded band is the power removed by trapping/destruction.
  Bottom row -- HOW it is done: the surviving fraction
                 P = (beta + P_el)(1 + y) / (beta + P_el + y + P_dest)
                 decomposed into its four competing channels -- line (Sobolev)
                 escape beta, electron-scattering escape P_el, continuum
                 destruction P_dest, and collisional quenching y = C_ul/A_ul.

Only continuum destruction (P_dest) is fed back as an explicit DAO heating term.
The collisional de-excitation rate y = C_ul/A_ul enters the survival
normalization but is not added as a separate DAO heating source.

Data (local to this folder; model 603b2ef4, latest iteration 025):
  line_escape_selected_603b2ef4.dat -- per-depth channel profiles for the four
    representative lines.
  line_escape_lines_603b2ef4.dat    -- per-line slab summary.
"""

import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "escape_probability.pdf")
LINES_FILE = os.path.join(HERE, "line_escape_lines_603b2ef4.dat")
SELECTED_FILE = os.path.join(HERE, "line_escape_selected_603b2ef4.dat")

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import LogLocator, NullFormatter

C_GRID = "#D6D9DE"


def selected_file_is_representative(path):
    with open(path, "r", encoding="utf-8") as fh:
        for _ in range(40):
            line = fh.readline()
            if not line:
                break
            if "representative escape-physics cases" in line:
                return True
    return False


def clean_line_label(label, energy_eV):
    label = " ".join(label.split())
    if not label:
        label = f"{energy_eV:.0f} eV"
    return f"{label} ({energy_eV:.0f} eV)"


def read_line_summary(path):
    rows = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 16:
                continue
            rows.append({
                "rank": int(parts[0]),
                "ip": int(parts[1]),
                "E_eV": float(parts[2]),
                "wave": float(parts[3]),
                "thin_sum": float(parts[4]),
                "escaped_sum": float(parts[5]),
                "heat_sum": float(parts[6]),
                "slab_ratio": float(parts[7]),
                "heat_ratio": float(parts[8]),
                "mean_beta": float(parts[9]),
                "mean_Pelec": float(parts[10]),
                "mean_Pdest": float(parts[11]),
                "mean_y": float(parts[12]),
                "mean_tauT": float(parts[13]),
                "max_tau_line": float(parts[14]),
                "label": parts[15],
                "comment": parts[16] if len(parts) > 16 else "",
            })
    return sorted(rows, key=lambda row: row["rank"])


def read_selected_profiles(path):
    groups = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 14:
                continue
            sel = int(parts[0])
            row = {
                "selected": sel,
                "depth": int(parts[1]),
                "tau": float(parts[2]),
                "ip": int(parts[3]),
                "E_eV": float(parts[4]),
                "thin": float(parts[5]),
                "escaped": float(parts[6]),
                "ratio": float(parts[7]),
                "heat": float(parts[8]),
                "beta": float(parts[9]),
                "Pelec": float(parts[10]),
                "Pdest": float(parts[11]),
                "y": float(parts[12]),
                "label": parts[13],
                "comment": parts[14] if len(parts) > 14 else "",
            }
            groups.setdefault(sel, []).append(row)

    profiles = []
    for sel in sorted(groups):
        rows = sorted(groups[sel], key=lambda r: r["tau"])
        if not rows:
            continue
        tau = np.array([r["tau"] for r in rows])
        ratio = np.array([r["ratio"] for r in rows])
        beta = np.array([r["beta"] for r in rows])
        pelec = np.array([r["Pelec"] for r in rows])
        pdest = np.array([r["Pdest"] for r in rows])
        yq = np.array([r["y"] for r in rows])
        thin = np.array([r["thin"] for r in rows])
        escaped = np.array([r["escaped"] for r in rows])
        E = rows[0]["E_eV"]
        label = clean_line_label(rows[0]["label"], E)
        thin_sum = np.sum(thin)
        escaped_sum = np.sum(escaped)
        profiles.append({
            "selected": sel,
            "ip": rows[0]["ip"],
            "E_eV": E,
            "raw_label": rows[0]["label"],
            "comment": rows[0].get("comment", ""),
            "tau": tau,
            "ratio": ratio,
            "beta": beta,
            "pelec": pelec,
            "pdest": pdest,
            "y": yq,
            "thin": thin,
            "escaped": escaped,
            "label": label,
            "thin_sum": thin_sum,
            "escaped_sum": escaped_sum,
            "zeros_thin": int(np.sum(thin == 0.0)),
            "zeros_escaped": int(np.sum(escaped == 0.0)),
            "slab_ratio": escaped_sum / thin_sum if thin_sum > 0 else 0.0,
        })
    return profiles


def selected_lines_from_diagnostic(lines_path, selected_path, n=4):
    summaries = {row["ip"]: row for row in read_line_summary(lines_path)}
    selected = []
    for profile in read_selected_profiles(selected_path)[:n]:
        line = dict(profile)
        summary = summaries.get(line["ip"])
        if summary is not None:
            line.update({
                "rank": summary["rank"],
                "label": clean_line_label(summary["label"], summary["E_eV"]),
                "thin_sum": summary["thin_sum"],
                "escaped_sum": summary["escaped_sum"],
                "slab_ratio": summary["slab_ratio"],
                "heat_ratio": summary["heat_ratio"],
                "max_tau_line": summary["max_tau_line"],
            })
        else:
            line.setdefault("rank", -1)
            line.setdefault("heat_ratio", 0.0)
            line.setdefault("max_tau_line", 0.0)
        selected.append(line)
    if len(selected) < n:
        raise RuntimeError(
            f"only found {len(selected)} selected line profiles in {selected_path}"
        )
    return selected


# ---------------------------------------------------------------------------
# Journal figure.
#
# Two-row small-multiple layout (no dual y-axes): the top row shows what the
# escape branching does to the emergent line power, the bottom row shows how it
# is done -- the survival fraction P decomposed into its four competing
# channels.  The four columns are the representative lines chosen by
# write_line_escape_diagnostics(), one per escape-physics regime.
# ---------------------------------------------------------------------------

# Semantic colours (colourblind-safe; from the validated data-viz palette).
C_ORIG = "#52514e"   # original / optically-thin line power (ink, dashed)
C_ESC = "#2a78d6"    # escaped line power  (blue)
C_REMOVED = "#8a8f98"  # line power not surviving as escaped emission
C_SURV = "#111417"   # survival fraction P (ink, heavy)
C_CH_BETA = "#2a78d6"   # beta  -- line (Sobolev) escape       (blue)
C_CH_ELEC = "#1baf7a"   # P_elec -- electron-scattering escape (aqua)
C_CH_DEST = "#e34948"   # P_dest -- continuum destruction      (red)
C_CH_Y = "#eb6834"      # y      -- collisional de-excitation  (orange)

_ROMAN = ["", "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
          "XI", "XII", "XIII", "XIV", "XV", "XVI", "XVII", "XVIII", "XIX",
          "XX", "XXI", "XXII", "XXIII", "XXIV", "XXV", "XXVI", "XXVII"]


def spectroscopic_label(raw_label, comment, E_eV):
    """Turn a Cloudy line label ('O  8   18.9689A') into 'O VIII 18.97 Å'."""
    m = re.match(r"\s*([A-Z][a-z]?)\s+(\d+)\s+([\d.]+)\s*([A-Za-z]*)",
                 raw_label or "")
    if not m:
        return clean_line_label(raw_label, E_eV)
    elem, stage = m.group(1), int(m.group(2))
    wave, unit = float(m.group(3)), m.group(4)
    roman = _ROMAN[stage] if 0 < stage < len(_ROMAN) else str(stage)
    unit = "Å" if unit.upper().startswith("A") else unit   # Angstrom
    lam = f"{wave:.2f}" if wave < 100 else f"{wave:.1f}"
    # H-/He-like Lyman-alpha resonance lines get the familiar spectroscopic tag.
    tag = ""
    if comment and ("1^2S -   2^2P" in comment or "1^1S -   2^1P" in comment):
        tag = " Lyα"   # Lyman-alpha
    return f"{elem} {roman}{tag}  {lam} {unit}"


plt.rcParams.update({
    "font.size": 9.5,
    "axes.linewidth": 0.8,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "mathtext.fontset": "cm",
    "axes.labelpad": 3.0,
})

lines_file = LINES_FILE
selected_file = SELECTED_FILE
is_representative_selected = selected_file_is_representative(selected_file)
selected_lines = selected_lines_from_diagnostic(lines_file, selected_file, n=4)

fig, axes = plt.subplots(
    2, 4, figsize=(11.0, 5.6), sharex=True, sharey="row",
    gridspec_kw=dict(height_ratios=[1.0, 1.0], hspace=0.10, wspace=0.06),
)
fig.subplots_adjust(left=0.08, right=0.985, top=0.865, bottom=0.135)
letters = ["(a)", "(b)", "(c)", "(d)"]
TINY = 1.0e-6   # probability floor for the log axis: below this a channel is off

for col, line in enumerate(selected_lines):
    tau = line["tau"]
    thin = line["thin"]
    esc = line["escaped"]

    # ------- top row: emergent line power -------------------------------
    axt = axes[0][col]
    mt = thin > 0.0
    me = esc > 0.0
    # Shade the power removed by the escape branching (thin -> escaped).
    both = mt & me
    axt.fill_between(tau[both], np.maximum(esc[both], TINY), thin[both],
                     where=thin[both] > esc[both], color=C_REMOVED,
                     alpha=0.22, lw=0, zorder=1, label="removed")
    axt.plot(tau[mt], thin[mt], color=C_ORIG, lw=1.6, ls=(0, (5, 2)),
             zorder=3, label="optically thin")
    axt.plot(tau[me], esc[me], color=C_ESC, lw=2.3, zorder=4, label="escaped")
    axt.set_yscale("log")
    axt.set_xscale("log")

    # Panel header: (letter) ion + wavelength ; energy + escape fraction.
    ion = spectroscopic_label(line.get("raw_label", ""),
                              line.get("comment", ""), line.get("E_eV", 0.0))
    axt.set_title(f"{letters[col]}  {ion}", loc="left", fontsize=10.5, pad=6)
    axt.text(0.05, 0.06,
             r"$%.0f\,$eV" % line.get("E_eV", 0.0) + "\n"
             + r"$\Sigma_{\rm esc}/\Sigma_{\rm thin}=%.2f$" % line["slab_ratio"],
             transform=axt.transAxes, fontsize=10.5, color=C_SURV,
             va="bottom", ha="left", linespacing=1.4)

    # ------- bottom row: survival fraction and its channels -------------
    # Only where the line actually emits (thin > 0) are the branching ratios
    # defined; elsewhere the diagnostic stores 0, so mask to avoid floor spikes.
    axb = axes[1][col]
    mp = thin > 0.0
    tp = tau[mp]
    axb.plot(tp, np.clip(line["ratio"][mp], TINY, 1.5), color=C_SURV, lw=2.4,
             zorder=6, label=r"$P$ (survival)")
    axb.plot(tp, np.clip(line["beta"][mp], TINY, None), color=C_CH_BETA, lw=1.5,
             zorder=5, label=r"$\beta$ (line escape)")
    axb.plot(tp, np.clip(line["pelec"][mp], TINY, None), color=C_CH_ELEC, lw=1.5,
             zorder=5, label=r"$P_{\rm el}$ ($e^-$ scattering)")
    axb.plot(tp, np.clip(line["pdest"][mp], TINY, None), color=C_CH_DEST, lw=1.5,
             zorder=5, label=r"$P_{\rm dest}$ (destruction)")
    axb.plot(tp, np.clip(line["y"][mp], TINY, None), color=C_CH_Y, lw=1.5,
             zorder=4, label=r"$y$ (coll. quench)")
    axb.axhline(1.0, color=C_GRID, lw=0.8, ls=(0, (1, 2)), zorder=1)
    axb.set_yscale("log")
    axb.set_xscale("log")
    axb.set_ylim(TINY, 2.0)

    for ax in (axt, axb):
        ax.grid(True, which="major", color=C_GRID, lw=0.6, zorder=0)
        ax.tick_params(which="both", top=True, right=True, labelsize=8.6)
        ax.xaxis.set_major_locator(LogLocator(numticks=6))
        ax.xaxis.set_minor_locator(
            LogLocator(subs=np.arange(2, 10) * 0.1, numticks=100))
        ax.xaxis.set_minor_formatter(NullFormatter())

# Full slab depth (tau_T up to ~5). These H-/He-like ions of O and Ne exist
# only in the photoionized tau_T < 0.05 skin; deeper the gas recombines to
# lower stages, so the lines have zero emissivity and the profiles are empty --
# it is the ions that stop, not the depth grid.
tau_lo = min(l["tau"][l["tau"] > 0].min() for l in selected_lines)
tau_hi = max(l["tau"].max() for l in selected_lines)
axes[0][0].set_xlim(tau_lo * 0.75, tau_hi * 1.1)   # shared via sharex

# Shared top-row power scale (labels only on column 1, via sharey="row").
allthin = np.concatenate([l["thin"][l["thin"] > 0] for l in selected_lines])
allesc = np.concatenate([l["escaped"][l["escaped"] > 0] for l in selected_lines])
axes[0][0].set_ylim(max(allesc.min() * 0.5, allthin.max() * 1e-6),
                    allthin.max() * 3.0)

axes[1][0].set_xlabel("")  # x label is set once, centred, below

for col in range(4):
    axes[1][col].set_xlabel(r"Thomson depth  $\tau_{\rm T}$", fontsize=10)
axes[0][0].set_ylabel(r"Local line power" "\n" r"[erg cm$^{-3}$ s$^{-1}$]",
                      fontsize=10)
axes[1][0].set_ylabel(r"Branching terms" "\n" r"in survival factor", fontsize=10)

# Two legends: power (top row) and channels (bottom row), as figure legends.
power_handles = [
    Line2D([0], [0], color=C_ORIG, lw=1.6, ls=(0, (5, 2)), label="optically thin"),
    Line2D([0], [0], color=C_ESC, lw=2.3, label="escaped"),
    Patch(facecolor=C_REMOVED, alpha=0.22, label="not surviving as line emission"),
]
fig.legend(handles=power_handles, ncol=3, frameon=False, fontsize=9,
           loc="upper center", bbox_to_anchor=(0.5, 0.975),
           handlelength=1.8, columnspacing=1.6)

chan_handles = [
    Line2D([0], [0], color=C_SURV, lw=2.4, label=r"$P$  survival fraction"),
    Line2D([0], [0], color=C_CH_BETA, lw=1.6, label=r"$\beta$  line escape"),
    Line2D([0], [0], color=C_CH_ELEC, lw=1.6, label=r"$P_{\rm el}$  $e^-$-scattering escape"),
    Line2D([0], [0], color=C_CH_DEST, lw=1.6, label=r"$P_{\rm dest}$  continuum destruction"),
    Line2D([0], [0], color=C_CH_Y, lw=1.6, label=r"$y$  collisional quench"),
]
fig.legend(handles=chan_handles, ncol=5, frameon=False, fontsize=9,
           loc="lower center", bbox_to_anchor=(0.5, 0.012),
           handlelength=1.8, columnspacing=1.5)

fig.savefig(OUT)
fig.savefig(OUT.replace(".pdf", ".png"), dpi=200)

print(f"wrote {OUT}")
print(f"lines {lines_file}")
print(f"selected {selected_file}")
if not is_representative_selected:
    print("WARNING: selected file was generated by the old strongest-line rule; "
          "rerun maindaocl to get representative profiles.")
for line in selected_lines:
    print("rank=%d ip=%d %s thin_sum=%.6e escaped_sum=%.6e escaped/thin=%.3f "
          "zeros(thin,escaped)=(%d,%d)" %
          (line["rank"], line["ip"], line["label"], line["thin_sum"],
           line["escaped_sum"], line["slab_ratio"], line["zeros_thin"],
           line["zeros_escaped"]))
