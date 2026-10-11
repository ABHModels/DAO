"""Plot the final DAO spectra against the archived PEXRAV reflection spectrum."""

import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
MU_VIEW = 0.7
E_LO, E_HI = 100.0, 1e6  # 0.1–1000 keV, in eV
ISOTROPIC = "11168045"
BEAMS = ("310ae399", "04b092df", "3b9a421c", "607c5c2f")
ANGLE_AVERAGED_KERNEL = "3c9a43af"  # mu_inc requested as 0.79

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 8,
    "axes.linewidth": 0.6,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "legend.fontsize": 6.2,
    "legend.frameon": False,
    "pdf.fonttype": 42,
})


def read_run(run_hash):
    with (HERE / f"params_{run_hash}.json").open() as stream:
        params = json.load(stream)
    path = HERE / f"emergent_{run_hash}.dat"
    mu = []
    with path.open() as stream:
        for line in stream:
            if not line.startswith("#"):
                break
            match = re.search(r"I\(mu=([-\d.]+)\)", line)
            if match:
                mu.append(float(match.group(1)))
    data = np.loadtxt(path, comments="#")
    mu = np.asarray(mu)
    if data.shape[1] != 3 + len(mu) or not np.all(np.diff(data[:, 0]) > 0):
        raise ValueError(f"Invalid emergent grid: {path}")
    if not np.all(np.isfinite(data)):
        raise ValueError(f"Non-finite emergent spectrum: {path}")
    positive = np.flatnonzero(mu > 0)
    mu_pos = mu[positive]
    if not (mu_pos[0] < MU_VIEW < mu_pos[-1]):
        raise ValueError(f"Viewing cosine {MU_VIEW} is outside outgoing grid: {path}")
    upper = np.searchsorted(mu_pos, MU_VIEW)
    lo, hi = positive[upper - 1], positive[upper]
    weight = (MU_VIEW - mu[lo]) / (mu[hi] - mu[lo])
    intensity = (1 - weight) * data[:, 3 + lo] + weight * data[:, 3 + hi]
    return params, data[:, 0], mu_pos, intensity


def band_normalize(energy, intensity):
    inside = (energy >= E_LO) & (energy <= E_HI)
    if np.count_nonzero(inside) < 2:
        raise ValueError("The normalization band is outside the energy grid")
    mean = np.trapezoid(intensity[inside], energy[inside]) / (E_HI - E_LO)
    if not np.isfinite(mean) or mean <= 0:
        raise ValueError("Invalid normalization mean")
    return intensity / mean


runs = {run_hash: read_run(run_hash)
        for run_hash in (ISOTROPIC, *BEAMS, ANGLE_AVERAGED_KERNEL)}
run_specific = {"hash", "time", "label", "incidence", "angsca"}
reference_params, energy, mu_nodes, _ = runs[ISOTROPIC]
for run_hash, (params, other_energy, other_mu, _) in runs.items():
    common_params = {key: value for key, value in params.items()
                     if key not in run_specific}
    common_reference = {key: value for key, value in reference_params.items()
                        if key not in run_specific}
    if common_params != common_reference:
        raise ValueError(f"Physical parameters differ in run {run_hash}")
    if not np.array_equal(other_energy, energy) or not np.array_equal(other_mu, mu_nodes):
        raise ValueError(f"Emergent grids differ in run {run_hash}")
if reference_params["incidence"] != -2 or not reference_params["angsca"]:
    raise ValueError("Unexpected isotropic run configuration")
for run_hash in BEAMS:
    if not runs[run_hash][0]["angsca"]:
        raise ValueError(f"Expected angle-dependent kernel for {run_hash}")
angle_averaged_params = runs[ANGLE_AVERAGED_KERNEL][0]
if angle_averaged_params["angsca"] or angle_averaged_params["incidence"] != 0.79:
    raise ValueError("Unexpected angle-averaged run configuration")

pexrav = np.loadtxt(HERE / "pexrav.dat")
pexrav_energy = pexrav[:, 0] * 1e3  # keV to eV
positive_ref = pexrav[:, 2] > 0
if not np.all(np.diff(pexrav_energy) > 0) or not np.any(positive_ref):
    raise ValueError("Invalid PEXRAV energy grid or reflection data")
ref_energy = pexrav_energy[positive_ref]
ref_intensity = pexrav[positive_ref, 2]
within_ref = (energy >= ref_energy[0]) & (energy <= ref_energy[-1])
pexrav_curve = np.full_like(energy, np.nan)
pexrav_curve[within_ref] = np.exp(np.interp(
    np.log(energy[within_ref]), np.log(ref_energy), np.log(ref_intensity)))
# Zero outside PEXRAV's positive support for the band integral; leave gaps in
# the plotted curve so no value is implied beyond the reference data.
pexrav_normalized = band_normalize(energy, np.nan_to_num(pexrav_curve, nan=0))
pexrav_normalized[~within_ref] = np.nan

fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(7.2, 2.8), sharey=True)
angle_averaged_mu = mu_nodes[np.argmin(abs(
    mu_nodes - angle_averaged_params["incidence"]))]
ax_a.loglog(energy, band_normalize(energy, runs[ANGLE_AVERAGED_KERNEL][3]),
            color="#D55E00", lw=0.6, ls="--",
            label=rf"DAO, avg. kernel ($\mu_{{\rm inc}}={angle_averaged_mu:.3f}$)")
ax_a.loglog(energy, pexrav_normalized, color="#444444", lw=0.6,
            label="PEXRAV")

colors = ("#0072B2", "#009E73", "#D55E00", "#7B3294")
ax_b.loglog(energy, band_normalize(energy, runs[ISOTROPIC][3]),
            color="#B2182B", lw=0.6, label="DAO, isotropic incidence")
for run_hash, color in zip(BEAMS, colors):
    params, _, _, intensity = runs[run_hash]
    requested_mu = params["incidence"]
    snapped_mu = mu_nodes[np.argmin(abs(mu_nodes - requested_mu))]
    ax_b.loglog(energy, band_normalize(energy, intensity), color=color,
                lw=0.6, label=rf"DAO, $\mu_{{\rm inc}}={snapped_mu:.3f}$")
    print(f"{run_hash}: requested mu_inc={requested_mu:.2f}, "
          f"nearest GL node={snapped_mu:.6f}")
ax_b.loglog(energy, pexrav_normalized, color="#444444", lw=0.6,
            label="PEXRAV, reflection only")

for letter, axis in zip("ab", (ax_a, ax_b)):
    axis.set_xlim(1e2, 1e6)
    axis.set_ylim(1e-3, 1e3)
    axis.set_xlabel("Energy (eV)")
    axis.text(0.035, 0.96, f"({letter})", transform=axis.transAxes,
              ha="left", va="top", fontsize=8, fontweight="bold")
ax_a.set_ylabel(r"$I_{\rm refl}/\langle I\rangle_{0.1-1000\,\rm keV}$")
ax_a.legend(loc="lower left")
ax_b.legend(loc="lower left", fontsize=5.5)
ax_a.text(0.97, 0.96,
          rf"$\Gamma={reference_params['Gamma']}$" + "\n"
          rf"$E_{{\rm cut}}={reference_params['E_cut']}\,\rm keV$" + "\n"
          rf"$\log n_{{\rm H}}={reference_params['nh']}$" + "\n"
          rf"$\log\xi={reference_params['zeta']}$" + "\n"
          rf"$\mu_{{\rm view}}={MU_VIEW}$",
          transform=ax_a.transAxes, ha="right", va="top", fontsize=6.3,
          bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", lw=0.4))
fig.subplots_adjust(left=0.10, right=0.98, bottom=0.17, top=0.97, wspace=0.08)
for extension, options in (("png", {"dpi": 600}), ("pdf", {})):
    destination = HERE / f"compare_pexrav_inc.{extension}"
    fig.savefig(destination, bbox_inches="tight", **options)
    print(f"Saved {destination}")
plt.close(fig)
