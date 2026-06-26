"""
compps benchmark generator + plotter (PyXspec), driven by a DAO run hash.

Workflow:
  1. Read results/<hash>/params.json (written by maindaocl).
  2. Read the emergent RT spectra results/<hash>/emergent_compps.dat. The RT
     solver emits intensities on a set of angle (mu) nodes; the OUTGOING nodes
     (mu>0) are the *inclination* cosines at which we observe the slab, all
     stored in this one file. In the compps test mode the solver uses the
     double-Gauss angle grid (init_angle_double_gauss), so with NA=10 these
     nodes coincide exactly with compPS's internal quadrature abscissas
     {0.0469, 0.2308, 0.5, 0.7692, 0.9531} -- requesting cosIncl at them makes
     compPS's angular interpolation exact (no interpolation error).
  3. Build the Xspec compPS (Poutanen & Svensson 1996) model for the SAME
     isothermal pure-scattering slab the RT test mode solved, taking
       kTe     <- kT_e        (electron temperature [keV])
       kTbb    <- kT_bb        (seed blackbody [keV])
       tau     <- tau_slab     (vertical Thomson optical depth)
       cosIncl <- inclination  (one compPS evaluation PER inclination node)
     with geom=1 (slab), cov_frac=1 (clean slab), rel_refl=0 (no reflection),
     Maxwellian electrons (Gmin,Gmax<1).

     NOTE 1: par["incidence"] is the corona ILLUMINATION (incidence) angle, NOT
     the observer inclination. It must NOT be used for cosIncl. cosIncl is the
     emergent-direction inclination, i.e. each outgoing mu node from the RT file.
     NOTE 2: compPS's cosIncl has a hard range [0.05, 0.95], so the two extreme
     double-Gauss nodes (mu=0.0469 and 0.9531) fall just outside and are clamped
     -- the only nodes where compPS's interpolation is not exact. Clamped nodes
     are flagged on stdout.
  4. Overlay compPS against the RT emergent output, one curve per inclination.

Figure (journal style), two panels:
  - Left:  2x2 sub-panels, one per inclination; in each, RT (solid) vs compPS
           (dashed) normalised by their own total flux -> shape-only comparison.
  - Right: limb-darkening law F(mu)/F(mu_max) vs mu, with F = mu*I the observed
           flux (the mu area-projection that compPS also applies).

Energy grid is 0.01-1000 keV with 1000 log bins to match the RT grid.

Self-contained: the fixed run's emergent_compps.dat + params.json sit next to
this script, and the figure is written into the same folder. No arguments.
Requires HEASoft/PyXspec for the compPS model:
  export HEADAS=/path/to/heasoft/arch && source $HEADAS/headas-init.sh
  python comppsGenerator.py

Author:      Yimin Huang
Affiliation: Fudan University
Email:       huangym23@m.fudan.edu.cn
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt
import xspec

# Self-contained: the RT data and the output figure live in this script's own
# folder, so paths resolve from any cwd. The run is fixed (no CLI argument);
# its emergent_compps.dat and params.json are distributed alongside the script.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = SCRIPT_DIR
RUN_HASH = "8c3e57d7"

PLOT_TYPE = "emodel"   # our model saves emodel (E * photons)
# compPS cosIncl hard range; nodes outside are clamped (interpolation not exact
# there). Note the 5-pt double-Gauss endpoints (0.0469, 0.9531) both fall just
# outside [0.05, 0.95], so only the middle three nodes are exactly representable.
COSINCL_MIN = 0.05
COSINCL_MAX = 0.95


# ------------------------------------------------------------------
# Save the loaded model spectrum to a QDP file via `wd`
# (mirrors DataComparison/DataGenerator/xspecplot.py).
# ------------------------------------------------------------------
def write_data(path, plot_type=PLOT_TYPE):
    # QDP's `wd` truncates long path arguments, so write to a short basename in
    # the cwd first and then move it to the (possibly long) destination path.
    tmp = "_compps_qdp_tmp.dat"
    for p in (path, tmp):
        if os.path.exists(p):
            os.remove(p)
    xspec.Plot.device = "/null"
    xspec.Plot.xAxis = "keV"
    xspec.Plot.add = True
    xspec.Plot.addCommand(f"wd {tmp}")
    xspec.Plot(plot_type)
    xspec.Plot.addCommand("exit")
    xspec.Plot.commands = ()
    if not os.path.exists(tmp):
        raise RuntimeError(f"QDP did not write {tmp}")
    os.replace(tmp, path)


def load_params():
    pjson = os.path.join(PROJECT_ROOT, "params.json")
    with open(pjson) as f:
        par = json.load(f)
    print(f"Loaded {pjson}")
    return par


def generate_compps(par, cos_incl, out_path):
    """Build compPS for a single observer inclination (cos_incl) and save its
    emodel spectrum. cos_incl is an emergent-direction cosine, NOT the corona
    incidence angle."""
    kTe  = float(par["kT_e"])
    kTbb = float(par["kT_bb"])
    tau  = float(par["tau_slab"])

    print(f"compps: kTe={kTe} keV  kTbb={kTbb} keV  tau={tau}  "
          f"cosIncl={cos_incl:.6f}")

    # Match the RT energy grid: 0.01-1000 keV, 1000 log bins.
    xspec.AllModels.setEnergies("0.01 1000 1000 log")
    xspec.AllModels.clear()
    m = xspec.Model("compPS")

    c = m.compPS
    c.kTe.values      = kTe       # 1  electron temperature [keV] (>1 for Maxwellian)
    c.EleIndex.values = 2.0       # 2  power-law index (dummy for Maxwellian)
    c.Gmin.values     = 0.1       # 3  <1 -> Maxwellian
    # c.Gmax.values     = 0.1       # 4  <1 -> Maxwellian
    c.kTbb.values     = kTbb      # 5  seed blackbody [keV]
    c.tau_y.values    = tau       # 6  >0 -> vertical Thomson optical depth
    c.geom.values     = 1.0       # 7  slab (bottom injection)
    c.HovR_cyl.values = 1.0       # 8  H/R (dummy for slab)
    # 9  cosine of OBSERVER inclination; clamp into compPS's hard range
    # [0.05, 0.95] (edge-on and face-on double-Gauss nodes fall just outside).
    cos_set = min(max(cos_incl, COSINCL_MIN), COSINCL_MAX)
    if cos_set != cos_incl:
        print(f"  NOTE: cosIncl={cos_incl:.6f} outside compPS hard range "
              f"[{COSINCL_MIN}, {COSINCL_MAX}]; clamped to {cos_set:.4f} "
              f"(interpolation NOT exact at this node)")
    c.cosIncl.values  = cos_set
    c.cov_frac.values = 1.0       # 10 clean slab
    c.rel_refl.values = 0.0       # 11 reflection off
    c.Fe_ab_re.values = 1.0       # 12 dummy
    c.Me_ab.values    = 1.0       # 13 dummy
    c.xi.values       = 0.0       # 14 dummy
    c.Tdisk.values    = 1.0e6     # 15 dummy
    c.Betor10.values  = -10.0     # 16 no relativistic blurring
    c.Rin.values      = 10.0      # 17 dummy
    c.Rout.values     = 1000.0    # 18 dummy
    c.Redshift.values = 0.0       # 19 redshift
    c.norm.values     = 1.0

    write_data(out_path, PLOT_TYPE)
    print(f"Saved compps spectrum: {out_path}")
    return cos_set   # the cosIncl actually used (may be clamped)


def read_qdp(path):
    """Read a QDP file, returning (E[keV], emodel) numeric columns 1 and 3."""
    E, Y = [], []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                e = float(parts[0]); y = float(parts[2])
            except ValueError:
                continue   # skip QDP header lines
            E.append(e); Y.append(y)
    return np.array(E), np.array(Y)


def read_rt_emergent(path):
    """Read results/<hash>/emergent_compps.dat.

    Returns (E[keV], incl_mus, F_by_mu) where incl_mus is the sorted array of
    OUTGOING (mu>0) inclination cosines and F_by_mu[k] is the OBSERVED flux at
    inclination incl_mus[k].

    Projection factor (verified against compps.f): the RT solver outputs the
    bare emergent intensity I(mu,E) at the surface, but compPS reports the flux
    an observer measures, which for a flat slab seen at inclination i includes
    the foreshortening of the projected emitting area (cos i = mu):
        F_obs(mu,E) = mu * I(mu,E)
    In compps.f the emergent spectrum is interpolated from SUIPL(mu)*UANG(mu)
    = I(mu)*mu (UANG is the cosine grid, AANG the quadrature weight). We
    therefore multiply by mu here so the comparison is on equal footing.
    In the per-inclination shape panels this energy-independent factor cancels
    under each curve's own normalisation; it only changes the limb-darkening
    law (without it the bare intensity is mildly limb-BRIGHTENED, with it the
    observed flux is limb-DARKENED, matching compPS)."""
    mus = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") and "I_emergent(mu=" in line:
                mus.append(float(line.split("mu=")[1].split(")")[0]))
    data = np.loadtxt(path, comments="#")
    E_eV = data[:, 0]
    # data columns: 0=E, 1=I_corona, 2=I_disk, 3..=I(mu_nm) for nm=0..NA-1
    mus = np.array(mus)
    pos = np.where(mus > 0)[0]
    order = pos[np.argsort(mus[pos])]      # ascending mu (edge-on -> face-on)
    incl_mus = mus[order]
    F_by_mu = [mu * data[:, 3 + nm] for mu, nm in zip(incl_mus, order)]
    print(f"RT: outgoing inclination cosines = "
          f"{', '.join(f'{m:.4f}' for m in incl_mus)}")
    return E_eV / 1.0e3, incl_mus, F_by_mu


def total_flux(E, Y):
    """Integrate a (positive) spectrum over energy for normalisation."""
    m = np.isfinite(E) & np.isfinite(Y) & (Y > 0)
    if m.sum() < 2:
        return 1.0
    return np.trapz(Y[m], E[m])


def plot_compare(run_hash, par, Er, incl_mus, I_by_mu, compps_by_mu,
                 clamped_mus=()):
    """Two-panel journal figure.

    Left:  2x2 sub-panels, one per inclination; RT vs compPS each normalised
           by its own total flux -> spectral-shape comparison.
    Right: limb-darkening law F(mu)/F(mu_max) vs mu (F = mu*I, observed flux).

    clamped_mus: inclination nodes whose compPS cosIncl was clamped into the
    hard range [0.05, 0.95] (interpolation not exact there); reported on stdout.
    """
    # Nature-style rcParams (mirror of benchmark/plot_compare_pexrav.py).
    plt.rcParams.update({
        'font.family':        'sans-serif',
        'font.sans-serif':    ['Helvetica', 'Arial', 'DejaVu Sans'],
        'font.size':           7,
        'axes.labelsize':      8,
        'axes.titlesize':      8,
        'axes.linewidth':      0.6,
        'xtick.labelsize':     7,
        'ytick.labelsize':     7,
        'xtick.direction':     'in',
        'ytick.direction':     'in',
        'xtick.top':           True,
        'ytick.right':         True,
        'xtick.major.size':    3.0,
        'ytick.major.size':    3.0,
        'xtick.minor.size':    1.8,
        'ytick.minor.size':    1.8,
        'xtick.major.width':   0.6,
        'ytick.major.width':   0.6,
        'legend.fontsize':     6.5,
        'legend.frameon':      False,
        'pdf.fonttype':        42,
        'ps.fonttype':         42,
    })

    n = len(incl_mus)
    cmap = plt.get_cmap("viridis")
    colors = [cmap(0.1 + 0.8 * i / max(n - 1, 1)) for i in range(n)]
    thetas = np.degrees(np.arccos(incl_mus))

    # The left 2x2 panel can only show four inclinations. With many outgoing
    # nodes (e.g. NA=20 -> 10 nodes with mu>0), pick four evenly spaced ones
    # spanning edge-on (smallest mu) to face-on (largest mu), always keeping
    # the two endpoints. The right-panel limb law still uses all n nodes.
    n_show = min(4, n)
    left_idx = np.unique(np.linspace(0, n - 1, n_show).round().astype(int))
    print(f"Left panel inclinations (mu): "
          f"{', '.join(f'{incl_mus[k]:.4f}' for k in left_idx)}")

    fig = plt.figure(figsize=(7.2, 3.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 0.85], wspace=0.30)
    gs_l = gs[0, 0].subgridspec(2, 2, hspace=0.0, wspace=0.0)
    axes_l = [fig.add_subplot(gs_l[i // 2, i % 2]) for i in range(4)]
    ax_right = fig.add_subplot(gs[0, 1])

    # consistent energy range for the left sub-panels
    e_lo = max(Er[Er > 0].min(), 1e-2)
    e_hi = Er.max()

    # ---- Left panel: 2x2 flux-normalised shape comparison per inclination ----
    for cell, idx in enumerate(left_idx):
        ax = axes_l[cell]
        mu = incl_mus[idx]
        Ir = I_by_mu[idx]
        Ec, Yc = compps_by_mu[mu]
        rt_norm = Ir / total_flux(Er, Ir)
        cp_norm = Yc / total_flux(Ec, Yc)
        ax.loglog(Er, rt_norm, lw=0.9, color=colors[idx], label="DAO")
        ax.loglog(Ec, cp_norm, ls="--", lw=0.9, color="k", alpha=0.8,
                  label="compPS")
        ax.set_xlim(e_lo, e_hi)
        ax.text(0.65, 0.93,
                rf"$\mu={mu:.3f}$" + "\n" + rf"$\theta={thetas[idx]:.0f}^\circ$",
                transform=ax.transAxes, fontsize=6.5, va="top", ha="left")
        ax.grid(True, which="both", alpha=0.2)
        row, col = cell // 2, cell % 2
        if row == 0:
            ax.tick_params(labelbottom=False)
        if col == 1:
            ax.tick_params(labelleft=False)
        if cell == 0:
            ax.legend(loc="lower left", borderpad=0.3)

    # Parameter box (pexrav style) in the empty lower-left sub-panel (cell 2).
    param_txt = (
        'compPS slab\n'
        rf'$kT_e = {par["kT_e"]}\,$keV' + '\n'
        rf'$kT_{{bb}} = {par["kT_bb"]}\,$keV' + '\n'
        rf'$\tau = {par["tau_slab"]}$'
    )
    axes_l[2].text(0.05, 0.40, param_txt, transform=axes_l[2].transAxes,
                   ha="left", va="top", fontsize=6.5,
                   bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7",
                             lw=0.4))

    axes_l[2].set_xlabel(r"$E\ \mathrm{[keV]}$")
    axes_l[3].set_xlabel(r"$E\ \mathrm{[keV]}$")
    fig.text(0.025, 0.55, r"$I(E)\,/\,\int I(E)\,dE$",
             va="center", rotation="vertical", fontsize=8)
    fig.text(0.30, 0.935, "Flux-normalised emergent spectra", ha="center",
             fontsize=8)

    # ---- Right panel: limb-darkening law F(mu)/F(mu_max) vs mu ----
    # Energy-integrated OBSERVED flux F = mu*I per inclination (the mu
    # projection is applied in read_rt_emergent so RT matches compPS, which
    # also reports mu*I). The ratio to the most face-on node is dimensionless,
    # so RT (physical units) and compPS (Xspec units) compare on the same axis.
    rt_int = np.array([total_flux(Er, I_by_mu[i]) for i in range(n)])
    cp_int = np.array([total_flux(*compps_by_mu[mu]) for mu in incl_mus])
    rt_law = rt_int / rt_int[-1]
    cp_law = cp_int / cp_int[-1]
    ax_right.plot(incl_mus, rt_law, "o-", lw=1.0, ms=4, color="C0",
                  label="DAO")
    ax_right.plot(incl_mus, cp_law, "s--", lw=0.9, ms=4, color="C3",
                  label="compPS (Xspec)")
    # x-ticks only at the actual angle nodes (mu values)
    ax_right.set_xticks(incl_mus)
    ax_right.set_xticklabels([f"{m:.3f}" for m in incl_mus])
    ax_right.set_xlabel(r"$\mu = \cos\,i$")
    ax_right.set_ylabel(r"$F(\mu)\,/\,F(\mu_\mathrm{max})$,"
                        r"  $F=\mu\,I$  (energy-integrated)")
    ax_right.legend(loc="upper left", borderpad=0.3)
    ax_right.grid(True, alpha=0.3)
    fig.text(0.79, 0.935, "Limb-darkening law", ha="center", fontsize=8)

    fig.subplots_adjust(left=0.10, right=0.97, bottom=0.22, top=0.91)

    # edge-on -> face-on directional cue: a clean straight arrow beneath the mu
    # axis spanning the right panel, decoupled from the data so it never
    # overlaps the curves. Placed in figure coordinates after layout is fixed.
    fig.canvas.draw()
    pos = ax_right.get_position()
    xL, xR, W = pos.x0, pos.x1, pos.width
    yrow = pos.y0 - 0.155              # just below the "mu = cos i" axis label
    TXT_C = "k"                        # text colour
    ARR_C = "0.35"                     # arrow colour (muted grey)
    ax_right.annotate("", xy=(xL + 0.66 * W, yrow), xytext=(xL + 0.27 * W, yrow),
                      xycoords="figure fraction", textcoords="figure fraction",
                      arrowprops=dict(arrowstyle="-|>", color=ARR_C, lw=1.0,
                                      mutation_scale=13))
    fig.text(xL + 0.01 * W, yrow, "edge-on", ha="left", va="center",
             fontsize=8, fontweight="bold", color=TXT_C)
    fig.text(xR - 0.01 * W, yrow, "face-on (centre)", ha="right", va="center",
             fontsize=8, fontweight="bold", color=TXT_C)

    plt.show()
    out_base = os.path.join(PROJECT_ROOT, f"compps_compare_{run_hash}")
    fig.savefig(f"{out_base}.png", dpi=200)
    fig.savefig(f"{out_base}.pdf", dpi=300)
    print(f"Saved plot: {out_base}.png and {out_base}.pdf")


def main():
    par = load_params()

    rt_path = os.path.join(PROJECT_ROOT, "emergent_compps.dat")
    Er, incl_mus, I_by_mu = read_rt_emergent(rt_path)

    # One compPS evaluation per inclination node (cosIncl = the RT node, which
    # with the double-Gauss grid is a compPS native abscissa). Track which nodes
    # had to be clamped to compPS's cosIncl hard max (requires HEASoft/PyXspec).
    compps_by_mu = {}
    clamped_mus = []
    for mu in incl_mus:
        out_path = os.path.join(PROJECT_ROOT,
                                f"compps_{RUN_HASH}_mu{mu:.4f}.dat")
        cos_used = generate_compps(par, mu, out_path)
        compps_by_mu[mu] = read_qdp(out_path)
        if abs(cos_used - mu) > 1e-9:
            clamped_mus.append(mu)

    plot_compare(RUN_HASH, par, Er, incl_mus, I_by_mu, compps_by_mu, clamped_mus)


if __name__ == "__main__":
    main()
