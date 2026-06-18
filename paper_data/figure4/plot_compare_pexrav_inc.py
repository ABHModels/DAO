import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import re
import json
import os

# ── Nature-style rcParams (mirror of plot_compare_pexrav.py) ─────────────────
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

# ── Configuration ────────────────────────────────────────────────────────────
# Self-contained: every input file sits next to this script, named by run hash:
#   emergent_<hash>.dat  (final-iteration emergent spectrum)
#   params_<hash>.json   (run parameters)
#   pexrav.dat           (pexrav reflection-only reference)
HERE = os.path.dirname(os.path.abspath(__file__))

# Panel (a): incidence-averaged spectrum, precomputed and stored as a run.
AVG_HASH = '38ea5924_inc_avg'

# Panel (a) extra: DAO v2.0 angle-averaged emergent intensity (last iteration).
DAO2_HASH = '4f606092'

# Panel (b): individual runs that differ only in the corona-illumination
# incidence angle μ_inc (cutoffpl, logξ=0, Γ=2, E_cut=300, nh=15).
RUN_HASHES = [
    '7743edbb',  # μ_inc = 0.1
    '9c6a048c',  # μ_inc = 0.2
    '5d9b4915',  # μ_inc = 0.3
    '1ecaa446',  # μ_inc = 0.4
    '514aadc7',  # μ_inc = 0.5
    '53250818',  # μ_inc = 0.6
    '439669b2',  # μ_inc = 0.8
    'b5e77f13',  # μ_inc = 0.9
]
MU_OBS_TARGET = 0.7   # observation/inclination angle of the emergent spectrum


# ── Helpers ──────────────────────────────────────────────────────────────────
def emergent_path(hash_):
    """Path to the (flattened) final-iteration emergent file for a run."""
    return os.path.join(HERE, f'emergent_{hash_}.dat')


def read_emergent(em_file):
    """Read emergent file → (E_eV, mu_arr, data)."""
    mu_vals = []
    with open(em_file) as f:
        for line in f:
            if not line.startswith('#'):
                break
            m = re.search(r'I\(mu=([-\d.]+)\)', line)
            if m:
                mu_vals.append(float(m.group(1)))
    mu_arr = np.array(mu_vals)
    data   = np.loadtxt(em_file, comments='#')
    return data[:, 0], mu_arr, data


def I_at_mu(mu_arr, data, mu_target):
    """Linear-in-μ interpolation of I(E, μ) at μ = mu_target over outgoing μ>0."""
    pos_idx = np.where(mu_arr > 0)[0]
    mu_pos  = mu_arr[pos_idx]
    j_hi    = np.searchsorted(mu_pos, mu_target)
    j_lo    = j_hi - 1
    lo_i, hi_i = pos_idx[j_lo], pos_idx[j_hi]
    mu_lo, mu_hi = mu_arr[lo_i], mu_arr[hi_i]
    w_hi = (mu_target - mu_lo) / (mu_hi - mu_lo)
    return (1 - w_hi) * data[:, 3 + lo_i] + w_hi * data[:, 3 + hi_i]


def I_angle_avg(mu_arr, data):
    """Angle-average of emergent I(E, μ) over the outgoing (μ>0) hemisphere.

    Uses Gauss-Legendre weights recovered from the run's μ nodes, so the result
    is the proper quadrature mean ∫₀¹ I dμ / ∫₀¹ dμ over the emergent directions.
    """
    nodes, weights = np.polynomial.legendre.leggauss(len(mu_arr))
    w = np.array([weights[np.argmin(np.abs(nodes - mu))] for mu in mu_arr])
    pos = np.where(mu_arr > 0)[0]
    num = sum(w[i] * data[:, 3 + i] for i in pos)
    den = sum(w[i] for i in pos)
    return num / den


def interp_log(E_target, E_src, F_src):
    return np.exp(np.interp(np.log(E_target),
                            np.log(E_src),
                            np.log(np.clip(F_src, 1e-300, None))))


# ── Load pexrav reflection (reference, shared by both panels) ────────────────
px = np.loadtxt(os.path.join(HERE, 'pexrav.dat'))
E_px_eV  = px[:, 0] * 1e3
F_px_ref = px[:, 2]                  # reflection only

# Normalisation band (same as the single-angle / avg scripts)
E_LO, E_HI = 0.1e3, 1000e3


def band_mean(I, E_eV):
    band = (E_eV >= E_LO) & (E_eV <= E_HI)
    return np.trapz(I[band], E_eV[band]) / (E_HI - E_LO)


# ── Panel (a) data: precomputed incidence average ───────────────────────────
avg_file = emergent_path(AVG_HASH)
E_eV, mu_arr, data_avg = read_emergent(avg_file)
I_avg = I_at_mu(mu_arr, data_avg, MU_OBS_TARGET)
with open(os.path.join(HERE, f'params_{AVG_HASH}.json')) as f:
    par = json.load(f)
print(f'Panel (a): incidence-averaged ← {os.path.basename(avg_file)}')

F_px_on_avg = interp_log(E_eV, E_px_eV, F_px_ref)
I_avg_n = I_avg       / band_mean(I_avg, E_eV)
F_px_n  = F_px_on_avg / band_mean(F_px_on_avg, E_eV)

# DAO v2.0: emergent intensity at the same observation angle as the blue curve
# (μ_obs = MU_OBS_TARGET). This run differs only in the Compton-scattering
# treatment, so the comparison isolates that effect at fixed inclination.
dao2_file = emergent_path(DAO2_HASH)
E_d2, mu_d2, data_d2 = read_emergent(dao2_file)
I_d2 = I_at_mu(mu_d2, data_d2, MU_OBS_TARGET)
if not (E_d2.shape == E_eV.shape and np.allclose(E_d2, E_eV)):
    I_d2 = interp_log(E_eV, E_d2, I_d2)
I_d2_n = I_d2 / band_mean(I_d2, E_eV)
print(f'Panel (a): DAO v2.0 (Compton scattering) ← {os.path.basename(dao2_file)}')
# ── Panel (b) data: reflected spectrum at each angle bin ────────────────────
# The code (snap_incidence, source/params.cpp) snaps the requested -incidence to
# the nearest GL node g.mu[], so the only physically distinct spectra are those
# at the outgoing GL nodes.  Group runs by their snapped node and plot one curve
# per real angle bin, labelled by the node value.
GL_NODES = np.sort(mu_arr[mu_arr > 0])   # outgoing GL nodes from the run header

bins = {}   # snapped-node index -> {'mu_bin', 'curve', 'matched', 'mus'}
for h in RUN_HASHES:
    em_file = emergent_path(h)
    E_h, mu_h, data_h = read_emergent(em_file)
    I_h = I_at_mu(mu_h, data_h, MU_OBS_TARGET)
    matched = (E_h.shape == E_eV.shape and np.allclose(E_h, E_eV))
    if not matched:
        I_h = interp_log(E_eV, E_h, I_h)
    with open(os.path.join(HERE, f'params_{h}.json')) as f:
        mu_req = json.load(f)['incidence']
    k = int(np.argmin(np.abs(mu_req - GL_NODES)))   # snapped bin (same rule as code)

    cur = bins.get(k)
    # Pick the representative run for this bin: prefer one on the standard energy
    # grid (no interpolation).
    if cur is None:
        bins[k] = {'mu_bin': GL_NODES[k], 'curve': I_h,
                   'matched': matched, 'mus': [mu_req]}
    else:
        cur['mus'].append(mu_req)
        if matched and not cur['matched']:
            cur.update(curve=I_h, matched=matched)

groups = [bins[k] for k in sorted(bins)]
print(f'{len(RUN_HASHES)} runs → {len(groups)} angle bins (GL nodes):')
for g in groups:
    flag = '' if g['matched'] else '  [WARNING: representative on non-standard grid]'
    print(f"  μ_bin={g['mu_bin']:.4f}  ← requested "
          f"{', '.join(f'{m:.3f}' for m in sorted(g['mus']))}{flag}")

# pexrav on the common grid, normalised the same way
F_px_b_n = F_px_on_avg / band_mean(F_px_on_avg, E_eV)

# ── Plot ─────────────────────────────────────────────────────────────────────
C_DAO = '#0072B2'   # blue
C_PEX = '#444444'   # dark grey
# Discrete, well-separated colors — one per distinct spectrum.
inc_colors = cm.turbo(np.linspace(0.05, 0.95, len(groups)))

fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(6.8, 2.7), sharey=True)

# Panel (a): incidence-averaged vs pexrav
ax_a.loglog(E_eV, I_avg_n, color=C_DAO, lw=0.8, alpha=0.9,
            label=r'DAOv2.0 (ang-dep Compton, inc. avg)', zorder=3)
ax_a.loglog(E_eV, I_d2_n, color='#D55E00', lw=0.8, alpha=0.9,
            label=r'DAOv2.0 (ang-avg Compton, '
                  r'$\mu_{\rm inc}=0.7$)', zorder=4)
ax_a.loglog(E_eV, F_px_n, color=C_PEX, lw=0.9, alpha=0.95,
            label='pexrav (reflection-only)', zorder=5)
ax_a.set_xlim(1e2, 1e6)
ax_a.set_ylim(1e-4, 1e5)
ax_a.set_xlabel(r'Energy (eV)')
ax_a.set_ylabel(r'$I_{\rm refl} / \langle I\rangle$')
ax_a.legend(loc='lower left', fontsize=6.5, borderpad=0.3)
ax_a.text(0.04, 0.96, '(a)', transform=ax_a.transAxes,
          ha='left', va='top', fontsize=8, fontweight='bold')

incident_txt = (
    'cutoffpl\n'
    rf'$\Gamma = {par["Gamma"]}$' + '\n'
    rf'$E_{{\rm cut}} = {par["E_cut"]}\,$keV' + '\n'
    rf'$n_{{\rm H}} = 10^{{{par["nh"]}}}\,$cm$^{{-3}}$' + '\n'
    rf'$\log\xi = {par["zeta"]}$' + '\n'
    r'$\mu_{\rm view} = \cos 45^\circ$'
)
ax_a.text(0.97, 0.96, incident_txt, transform=ax_a.transAxes,
          ha='right', va='top', fontsize=6.5,
          bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='0.7', lw=0.4))

# Panel (b): reflected spectrum at each GL angle bin vs pexrav
for g, c in zip(groups, inc_colors):
    I_n = g['curve'] / band_mean(g['curve'], E_eV)
    lbl = rf'$\mu_{{\rm inc}}={g["mu_bin"]:.3f}$'
    ax_b.loglog(E_eV, I_n, color=c, lw=0.9, alpha=0.95, label=lbl, zorder=3)
ax_b.loglog(E_eV, F_px_b_n, color=C_PEX, lw=1.2, alpha=0.95,
            label='pexrav (refl.)', zorder=5)
ax_b.set_xlim(1e2, 1e6)
ax_b.set_xlabel(r'Energy (eV)')
ax_b.legend(loc='lower left', fontsize=5.8, borderpad=0.3,
            ncol=2, columnspacing=1.0, labelspacing=0.25)
ax_b.text(0.04, 0.96, '(b)', transform=ax_b.transAxes,
          ha='left', va='top', fontsize=8, fontweight='bold')

plt.subplots_adjust(left=0.09, right=0.98, bottom=0.15, top=0.97, wspace=0.08)

plt.show()
out_png = os.path.join(HERE, 'compare_pexrav_inc.png')
out_pdf = os.path.join(HERE, 'compare_pexrav_inc.pdf')
fig.savefig(out_png, dpi=600, bbox_inches='tight')
fig.savefig(out_pdf,            bbox_inches='tight')
print(f'Saved: {out_png}')
print(f'Saved: {out_pdf}')
