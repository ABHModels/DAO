"""Reproduce the all-line P=1, thermally re-equilibrated spectral comparison.
Run after collect.py has copied this script and the converged data to paper_data.
The two DAO curves share normalization factors to expose their actual change.
"""
from pathlib import Path
import hashlib
import json
import os
import tempfile

os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'dao-oxygen-highdensity-mpl'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, LogLocator, NullFormatter

HERE = Path(__file__).resolve().parent
OUT = HERE
(OUT / 'output/pdf').mkdir(parents=True, exist_ok=True)
(OUT / 'data').mkdir(exist_ok=True)
INPUTS = {
    'DAO': (HERE / 'input/DAO_standard.dat', 0, 1),
    'DAO (all P=1)': (HERE / 'input/DAO_all_P1.dat', 0, 1),
    'reflionx': (HERE / 'input/reflionx_flux_cgs.dat', 0, 1),
    'xillverCp': (HERE / 'input/xillverCp_angle_average.dat', 2, 4),
}
COLORS = {'DAO': '#0072B2', 'DAO (all P=1)': '#CC79A7', 'reflionx': '#444444', 'xillverCp': '#D55E00'}
STYLES = {'DAO': '-', 'DAO (all P=1)': '--', 'reflionx': '--', 'xillverCp': '-.'}
spectra = {}
for name, (path, ec, fc) in INPUTS.items():
    a = np.loadtxt(path)
    e, f = a[:, ec], a[:, fc]
    assert np.all(np.isfinite(e)) and np.all(np.diff(e) > 0)
    assert np.all(np.isfinite(f)) and np.all(f >= 0)
    spectra[name] = (e, f)


def log_interp(target, e, f):
    return np.exp(np.interp(np.log(target), np.log(e),
                            np.log(np.clip(f, 1e-300, None))))


plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 9.,
    'axes.labelsize': 10., 'axes.linewidth': .8,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.top': True, 'ytick.right': True,
    'xtick.major.size': 4., 'ytick.major.size': 4.,
    'xtick.minor.size': 2., 'ytick.minor.size': 2.,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    'savefig.facecolor': 'white',
})
fig, (left, right) = plt.subplots(1, 2, figsize=(9.2, 3.75))
fig.subplots_adjust(left=.077, right=.985, bottom=.17, top=.82, wspace=.30)
inset = right.inset_axes([.51, .53, .455, .40])
inset.set_facecolor('white')
left.axvspan(20., 50., color='.91', lw=0, zorder=0)
grid = spectra['DAO'][0]
band = (grid >= 20.) & (grid <= 50.)
records = {}
shared_divisors = {}
for name, (e, f) in spectra.items():
    on_grid = f if name == 'DAO' else log_interp(grid, e, f)
    efe = e * f
    broad_divisor = np.trapezoid(grid[band] * on_grid[band], grid[band]) / 30.
    if name == 'DAO (all P=1)':
        broad_divisor = shared_divisors['broad']
    broad = efe / broad_divisor
    soft_divisor = 2.25 * log_interp(2.25, e, f)
    if name == 'DAO (all P=1)':
        soft_divisor = shared_divisors['soft']
    soft = e * f / soft_divisor
    inset_divisor = np.interp([.55, .75], e, efe).mean()
    if name == 'DAO':
        shared_divisors = dict(broad=broad_divisor, soft=soft_divisor, inset=inset_divisor)
    if name == 'DAO (all P=1)':
        inset_divisor = shared_divisors['inset']
    sf = efe / inset_divisor
    displayed = (e >= .525) & (e <= .8)
    if name != 'DAO (all P=1)':
        assert np.isclose(np.interp([.55, .75], e, sf).mean(), 1., rtol=1e-12)
    style = dict(color=COLORS[name], ls=STYLES[name], lw=1.05)
    left.loglog(e, np.ma.masked_less_equal(broad, 0.), label=(r'DAO (all $P=1$)' if name == 'DAO (all P=1)' else name), **style)
    right.loglog(e, np.ma.masked_less_equal(soft, 0.), **style)
    inset.loglog(e, np.ma.masked_less_equal(sf, 0.), **style)
    path = INPUTS[name][0]
    records[name] = {
        'input': str(path.relative_to(HERE)),
        'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        'broadband_EF_E_divisor': float(broad_divisor),
        'soft_EF_E_divisor': float(soft_divisor),
        'inset_native_EF_E_divisor': float(inset_divisor),
    }
    np.savetxt(OUT / 'data' / f'native_{name}.dat',
               np.column_stack((e, broad, soft)), fmt='%.16e',
               header='E_keV EF_E_over_20_50_keV_mean EF_E_over_EF_E_at_2p25_keV')
    np.savetxt(OUT / 'data' / f'inset_{name}.dat',
               np.column_stack((e[displayed], sf[displayed])), fmt='%.16e',
               header='E_keV native_EF_E_over_mean_at_0p55_and_0p75_keV; '
                      'native sampling; no smoothing or continuum subtraction')

left.set(xlim=(.1, 1000.), ylim=(1e-2, 1e2), xlabel='Energy (keV)',
         ylabel=r'$E F_E\,/\,\langle E F_E\rangle_{20\mathrm{-}50\,{\rm keV}}$')
right.set(xlim=(.3, 2.3), ylim=(.8, 40.), xlabel='Energy (keV)',
          ylabel=r'$E F_E\,/\,(E F_E)_{2.25\,{\rm keV}}$')
inset.set(xlim=(.525, .8), ylim=(.5, 12.))
for ax, xticks, yticks in ((right, [.5, 1., 2.], [1., 10.]),
                            (inset, [.6, .7], [1., 5., 10.])):
    ax.xaxis.set_major_locator(FixedLocator(xticks))
    ax.yaxis.set_major_locator(FixedLocator(yticks))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:g}'))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:g}'))
    ax.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)))
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_minor_formatter(NullFormatter())
inset.tick_params(labelsize=7., pad=2., which='both', length=2.5)
inset.set_ylabel(r'$\widehat{E F_E}$', fontsize=8., labelpad=1.)
for ax, label in ((left, '(a) Full spectrum'), (right, '(b) Soft band')):
    ax.text(0., 1.04, label, transform=ax.transAxes,
            ha='left', va='bottom', fontsize=9.2)
left.text(.04, .94, r'$\log\xi=3,\quad n_{\rm H}=10^{18}\,{\rm cm}^{-3}$'
          + '\n' + r'$\Gamma=1.8,\quad kT_{\rm e}=60\,{\rm keV}$',
          transform=left.transAxes, fontsize=8.4, va='top', linespacing=1.55)
handles, labels = left.get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(.54, 1.01),
           ncol=4, frameon=False, handlelength=2.5, columnspacing=1.4, fontsize=9.5)
stem = 'all_survival_off_comparison'
fig.savefig(OUT / f'{stem}.png', dpi=500)
fig.savefig(OUT / f'{stem}.svg')
fig.savefig(OUT / 'output/pdf' / f'{stem}.pdf')
plt.close(fig)

manifest = {
    'inputs': records,
    'reference': 'https://academic.oup.com/mnras/article/543/3/2633/8268889',
    'reference_figure': 10,
    'left': {'xlim_keV': [.1, 1000.], 'ylim': [1e-2, 1e2],
             'normalization': 'Mean E F_E over 20--50 keV, integrated on the DAO grid'},
    'right': {'xlim_keV': [.3, 2.3], 'ylim': [.8, 40.],
              'ordinate': 'E F_E / [E F_E at 2.25 keV]',
              'normalization': 'Reference models and standard DAO aligned at 2.25 keV; P=1 DAO uses the standard DAO divisor'},
    'inset': {'xlim_keV': [.525, .8], 'ylim': [.5, 12.],
              'normalization': 'Native E F_E divided by mean at 0.55 and 0.75 keV; both DAO curves share the standard DAO divisor',
              'endpoint_interpolation': 'Linear interpolation to the two normalization energies only',
              'sampling': 'Native energy samples, connected by straight segments',
              'smoothing': False},
    'all_axes': 'Logarithmic energy and flux',
    'main_sampling': 'Native grids, no smoothing',
    'model_runs': 'DAO standard 5b7fe321; DAO all P=1 fe48cbde, with full thermal convergence; reference models from previous comparison',
    'control_normalization': 'Both DAO curves share the standard DAO divisors in every panel',
    'continuum_subtraction': False, 'equivalent_width': False,
}
(OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
(OUT / 'caption.tex').write_text(r"""\caption{Effect of removing the line-survival correction at
$\log\xi=3$, $n_{\rm H}=10^{18}\,\mathrm{cm}^{-3}$, and $\Gamma=1.8$.
The standard \textsc{DAO} model is compared with a fully recomputed atmosphere
with $P_{u\ell}=1$ for all lines, alongside \textsc{reflionx} and \textsc{xillverCp}.
Left: $E F_E$ normalized by the mean over 20--50\,keV.
Right: soft-band spectra normalized at 2.25\,keV; the inset shows
O\,\textsc{viii}, normalized by the mean at 0.55 and 0.75\,keV.
In every panel, both \textsc{DAO} curves share the standard \textsc{DAO}
normalization. No smoothing is applied.}
\label{fig:all_survival_off_comparison}
""")
print(json.dumps({'figure': str(OUT / f'{stem}.png'), 'inputs': records}, indent=2))
