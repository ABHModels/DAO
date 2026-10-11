"""Compare saved local Fe K-alpha production at log xi = 1 and 2."""
from pathlib import Path
import hashlib
import json
import os
import sys
import tempfile

os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'dao-kalpha-paper-mpl'))
import numpy as np
import matplotlib as mpl
if '--show' not in sys.argv:
    mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, NullFormatter, MultipleLocator

HERE = Path(__file__).resolve().parent
COLORS = {'all': '#252525', 'FeI_V': '#8E44AD', 'FeXVII_XXII': '#008F83'}
LABELS = {'all': r'All Fe K$\alpha$', 'FeI_V': 'Fe I-V', 'FeXVII_XXII': 'Fe XVII-XXII'}
STYLES = {1: '-', 2: (0, (4.0, 2.3))}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_run(record):
    for name, entry in record['files'].items():
        assert sha(HERE / name) == entry['sha256'], name
    base = HERE / 'data' / f"logxi{record['log_xi']}"
    a = np.loadtxt(base / 'depth_power.dat')
    r = np.loadtxt(base / 'FeXVII_XXII_Kalpha_records.dat')
    summary = json.loads((base / 'analysis.json').read_text())
    assert a.shape == (48, 19) and np.isfinite(a).all()
    assert r.shape[1] == 9 and np.isfinite(r).all()
    assert np.all((r[:, 1] >= 17) & (r[:, 1] <= 22))
    assert np.all(r[:, 2] == 0) and np.all(np.isin(r[:, 3], [1, 2]))
    assert np.all(r[:, 0] == r[:, 0].astype(int))
    assert np.all((r[:, 0] >= 0) & (r[:, 0] < len(a)))
    assert np.all(r[:, 8] >= 0)
    tau, dz = a[:, 1], a[:, 2]
    assert np.all(np.diff(tau) > 0) and np.all(dz > 0)
    curves = {
        'all': a[:, 13], 'FeI_V': a[:, 10],
        'FeXVII_XXII': np.bincount(r[:, 0].astype(int), weights=r[:, 8], minlength=len(a)),
    }
    assert np.all(curves['all'] > 0) and np.all(curves['FeI_V'] >= 0)
    assert np.all(curves['FeI_V'] + curves['FeXVII_XXII'] <= curves['all'] * (1 + 1e-12))
    np.testing.assert_allclose(curves['FeI_V'] / curves['all'], a[:, 14], rtol=1e-13)
    powers = {name: float(eps @ dz) for name, eps in curves.items()}
    np.testing.assert_allclose(powers['all'], summary['all_Fe_Kalpha_production_per_area'], rtol=1e-13)
    np.testing.assert_allclose(powers['FeI_V'], summary['FeI_V_Kalpha_production_per_area'], rtol=1e-13)
    np.testing.assert_allclose(powers['FeXVII_XXII'] / powers['all'],
                               summary['fluorescence_parent_groups_fraction_of_all_Kalpha']['Fe_17_22'], rtol=1e-13)
    values = {
        'log_xi': record['log_xi'], 'source_run': record['source_run'],
        'accepted_iteration': record['accepted_iteration'],
        'power_per_area_erg_cm2_s': powers,
        'slab_fraction_percent': {name: 100 * powers[name] / powers['all'] for name in ['FeI_V', 'FeXVII_XXII']},
    }
    return tau, curves, values


def main():
    snapshot = json.loads((HERE / 'source_snapshot.json').read_text())
    loaded = {r['log_xi']: load_run(r) for r in snapshot['runs']}
    mpl.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 9,
        'axes.labelsize': 10, 'axes.linewidth': .8,
        'xtick.labelsize': 9, 'ytick.labelsize': 9,
        'xtick.direction': 'in', 'ytick.direction': 'in',
        'xtick.top': True, 'ytick.right': True, 'legend.frameon': False,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'savefig.facecolor': 'white',
    })
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(3.5, 4.35), sharex=True,
                                     gridspec_kw={'height_ratios': [1.65, 1]})
    fig.subplots_adjust(left=.205, right=.975, bottom=.14, top=.975, hspace=.075)
    for xi in [1, 2]:
        tau, curves, _ = loaded[xi]
        for name, eps in curves.items():
            # Mask true zeros only for display on a logarithmic ordinate.
            top.plot(tau, np.ma.masked_less_equal(eps, 0), color=COLORS[name],
                     ls=STYLES[xi], lw=1.55, zorder=2 if name == 'all' else 3)
            if name != 'all':
                bottom.plot(tau, 100 * eps / curves['all'], color=COLORS[name],
                            ls=STYLES[xi], lw=1.55)
    top.set_yscale('log')
    top.set_ylim(1e-4, 3e5)
    top.set_ylabel(r'$\epsilon_{\mathrm{K}\alpha}$ (erg cm$^{-3}$ s$^{-1}$)', labelpad=5)
    top.yaxis.set_major_locator(LogLocator(base=10, numticks=6))
    top.yaxis.set_minor_locator(LogLocator(base=10, subs=(2, 5), numticks=50))
    top.yaxis.set_minor_formatter(NullFormatter())
    # Separate legends make color (ion group) and line style (log xi) independent.
    top.legend(handles=[Line2D([], [], color=COLORS[name], lw=1.6, label=LABELS[name])
                        for name in COLORS],
               loc='upper left', bbox_to_anchor=(.035, .81), fontsize=9,
               borderaxespad=0, handlelength=2, handletextpad=.65, labelspacing=.5)
    bottom.legend(handles=[Line2D([], [], color='#333333', lw=1.6, ls=STYLES[xi],
                                  label=rf'$\log\xi={xi}$') for xi in [1, 2]],
                  loc='center left', bbox_to_anchor=(.035, .48), fontsize=9,
                  borderaxespad=0, handlelength=2.6, handletextpad=.65, labelspacing=.6)
    bottom.set_ylim(-3, 105)
    bottom.set_ylabel('Fraction of total (%)', labelpad=5)
    bottom.set_xlabel(r'$\tau_{\rm T}$', labelpad=5)
    bottom.yaxis.set_major_locator(MultipleLocator(25))
    bottom.yaxis.set_minor_locator(MultipleLocator(5))
    for ax in (top, bottom):
        ax.set_xscale('log')
        ax.set_xlim(min(v[0].min() for v in loaded.values()), 5)
        ax.xaxis.set_major_locator(LogLocator(base=10, numticks=6))
        ax.xaxis.set_minor_locator(LogLocator(base=10, subs=(2, 5), numticks=40))
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.tick_params(which='major', length=3.8, width=.75, pad=3)
        ax.tick_params(which='minor', length=2, width=.6)
    out = HERE / 'output'
    (out / 'pdf').mkdir(parents=True, exist_ok=True)
    fig.savefig(out / 'pdf/local_Kalpha_production_logxi12.pdf')
    fig.savefig(out / 'local_Kalpha_production_logxi12.png', dpi=600)
    values = {
        'runs': [loaded[xi][2] for xi in [1, 2]],
        'line_styles': {'log_xi_1': 'solid', 'log_xi_2': 'dashed'},
        'colors': COLORS,
        'upper_panel': 'Absolute local all-direction K-alpha production, before transfer; no scaling or smoothing',
        'lower_panel': '100 * epsilon_group / epsilon_all at the same depth and log xi',
        'slab_fraction': 'sum(epsilon_group * dz) / sum(epsilon_all * dz)',
        'zero_handling': 'True zeros masked only in the logarithmic upper panel',
        'figure_size_inches': [3.5, 4.35], 'png_dpi': 600,
    }
    (HERE / 'plot_values.json').write_text(json.dumps(values, indent=2) + '\n')
    artifacts = [p for p in HERE.rglob('*') if p.is_file() and p.name != 'artifact_hashes.json']
    (HERE / 'artifact_hashes.json').write_text(json.dumps(
        {str(p.relative_to(HERE)): sha(p) for p in sorted(artifacts)}, indent=2) + '\n')
    print(json.dumps(values['runs'], indent=2))
    print(out / 'pdf/local_Kalpha_production_logxi12.pdf')
    if '--show' in sys.argv:
        plt.show()
    plt.close(fig)


if __name__ == '__main__':
    main()
