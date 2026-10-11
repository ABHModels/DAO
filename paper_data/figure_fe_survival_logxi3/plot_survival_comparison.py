"""Single-column Fe XXV comparison; plot saved, validated transfer results only."""
from pathlib import Path
import hashlib
import json
import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/dao-fe-control-matplotlib")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"


def main():
    assert json.loads((DATA / "completed.json").read_text())["state"] == "completed"
    original = np.loadtxt(DATA / "full_survival.dat")
    changed = np.loadtxt(DATA / "full_wx_P1.dat")
    xillver = np.loadtxt(DATA / "logxi3_xillvercp.dat")
    assert np.array_equal(original[:, 0], changed[:, 0])
    for a in (original, changed, xillver):
        assert np.isfinite(a).all()
        assert np.all(np.diff(a[:, 0]) > 0)
    normalization = json.loads((DATA / "normalization.json").read_text())
    divisor = normalization["DAO_normalization_divisor_shared_by_both_curves"]
    energy = original[:, 0]
    dao = original[:, 1] / divisor
    dao_p1 = changed[:, 1] / divisor
    peak_mask = (energy >= 6.2) & (energy <= 7.1)
    peak_increase = 100 * (dao_p1[peak_mask].max() / dao[peak_mask].max() - 1)
    report = json.loads((DATA / "control_analysis.json").read_text())
    assert abs(peak_increase - report["peak"]["total_peak_increase_percent"]) < 1e-9

    output = HERE / "output"
    (output / "pdf").mkdir(parents=True, exist_ok=True)
    with plt.rc_context({
        "font.family": "DejaVu Sans", "font.size": 9.5,
        "axes.labelsize": 11, "axes.linewidth": .8,
        "xtick.labelsize": 9, "ytick.labelsize": 9,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "legend.frameon": False, "pdf.fonttype": 42, "ps.fonttype": 42,
        "savefig.facecolor": "white",
    }):
        fig, ax = plt.subplots(figsize=(3.5, 3.1))
        fig.subplots_adjust(left=.155, right=.98, bottom=.16, top=.97)
        ax.plot(xillver[:, 0], xillver[:, 1], color="#D55E00", lw=1.5,
                ls=(0, (4, 2)), label="XILLVER", zorder=2)
        ax.plot(energy, dao, color="#0072B2", lw=1.65, label="DAO", zorder=3)
        ax.plot(energy, dao_p1, color="#333333", lw=1.35,
                ls=(0, (1.5, 1.5)), label=r"DAO ($P=1$)", zorder=4)
        ax.set(xlim=(6.2, 7.1), ylim=(0, 7.65),
               xlabel="Energy (keV)", ylabel=r"$\widehat{F}_E$")
        ax.xaxis.set_major_locator(MultipleLocator(.2))
        ax.xaxis.set_minor_locator(AutoMinorLocator(4))
        ax.yaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.tick_params(which="major", length=4, width=.75)
        ax.tick_params(which="minor", length=2, width=.6)
        ax.legend(loc="upper left", fontsize=9, handlelength=2.3,
                  labelspacing=.45, borderaxespad=.6)
        ax.text(.965, .97, r"$\log\xi=3$", transform=ax.transAxes,
                ha="right", va="top", fontsize=10)
        assert len(ax.lines) == 3
        fig.savefig(output / "pdf/fe_xxv_model_comparison.pdf")
        fig.savefig(output / "fe_xxv_model_comparison.png", dpi=350)
        plt.close(fig)

    inputs = [DATA / "full_survival.dat", DATA / "full_wx_P1.dat",
              DATA / "logxi3_xillvercp.dat", DATA / "normalization.json"]
    manifest = {
        "curves": ["XILLVER", "DAO", "DAO (P=1)"],
        "P1_label_definition": "P_w=P_x=1 only; all other lines retain their original survival factors.",
        "panels": 1, "size_inches": [3.5, 3.1], "log_xi": 3,
        "sampling": "Native energy samples; no smoothing or rebinning.",
        "ordinate": "Full upper-face flux, with the existing 5/8 keV endpoint normalization; no continuum subtraction.",
        "DAO_normalization_divisor_shared_by_both_curves": divisor,
        "XILLVER_source": "Existing paper XILLVER-CP reference, retaining its recorded normalization.",
        "control": "In every depth cell, add epsilon_w*(1-P_w)/(4*pi*dE) and epsilon_x*(1-P_x)/(4*pi*dE) to the original material source, then solve full-slab transfer with the fixed converged atmosphere and original illumination.",
        "thermal_structure_recomputed": False,
        "full_spectral_peak_increase_percent": peak_increase,
        "input_sha256": {str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
    }
    (output / "fe_xxv_model_comparison_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"curves": manifest["curves"], "peak_increase_percent": peak_increase}))


if __name__ == "__main__":
    main()
