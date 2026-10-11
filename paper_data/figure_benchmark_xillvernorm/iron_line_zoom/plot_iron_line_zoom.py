"""Fe-K zoom: complete endpoint-normalized spectra with inset ratios to DAO."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile

os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'dao-iron-matplotlib'))
import numpy as np
import matplotlib as mpl
if '--show' not in sys.argv:
    mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MaxNLocator, ScalarFormatter

HERE = Path(__file__).resolve().parent
INPUTS = HERE.parent
NAMES = ['DAO', 'reflionx', 'xillvercp']
COLORS = {'DAO': '#0072B2', 'reflionx': '#444444', 'xillvercp': '#D55E00'}
STYLES = {'DAO': '-', 'reflionx': '--', 'xillvercp': '-.'}
PAPER = 'https://arxiv.org/html/2509.13411v1#S3.F10'
TRAPEZOID = getattr(np, 'trapezoid', None) or np.trapz


def interp(e, flux, energies):
    # Local linear interpolation preserves integrated narrow-line power under rebinning.
    assert np.all(np.diff(e) > 0)
    assert np.min(energies) >= e[0] and np.max(energies) <= e[-1]
    return np.interp(energies, e, flux)


def continuum(energies, anchors, fluxes):
    slope = np.log(fluxes[1] / fluxes[0]) / np.log(anchors[1] / anchors[0])
    return fluxes[0] * (np.asarray(energies) / anchors[0]) ** slope


def bin_average(e, f, edges):
    values = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        x = np.r_[lo, e[(e > lo) & (e < hi)], hi]
        y = interp(e, f, x)
        values.append(TRAPEZOID(y, x) / (hi - lo))
    return np.asarray(values)


def continuum_average(edges, anchors, fluxes):
    slope = np.log(fluxes[1] / fluxes[0]) / np.log(anchors[1] / anchors[0])
    amplitude = fluxes[0] / anchors[0] ** slope
    lo, hi = edges[:-1], edges[1:]
    if abs(slope + 1) < 1e-10:
        return amplitude * np.log(hi / lo) / (hi - lo)
    return amplitude * (hi ** (slope + 1) - lo ** (slope + 1)) / ((slope + 1) * (hi - lo))


def load_models(record):
    d = np.loadtxt(INPUTS / record['DAO_file'])
    mu, wt = np.polynomial.legendre.leggauss(8)
    assert d.shape[1] == 11
    f = 2 * np.pi * np.sum(d[:, 3:][:, mu > 0] * (mu * wt)[mu > 0], axis=1) * 1000
    r = np.loadtxt(INPUTS / record['reflionx_file'])
    x = np.loadtxt(INPUTS / record['xillver_file'])
    return {'DAO': (d[:, 0] / 1000, f),
            'reflionx': (r[:, 0], r[:, 1] / r[:, 0]),
            'xillvercp': (x[:, 0], x[:, 2] / x[:, 0])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--anchors', nargs=2, type=float, default=[5., 8.], metavar=('LOW', 'HIGH'))
    parser.add_argument('--show', action='store_true')
    args = parser.parse_args()
    anchors = np.asarray(args.anchors)
    assert 0 < anchors[0] < anchors[1]
    output = HERE / 'output'
    (output / 'pdf').mkdir(parents=True, exist_ok=True)
    (HERE / 'data').mkdir(exist_ok=True)
    records = sorted(json.loads((INPUTS / 'benchmark_inputs.json').read_text())['runs'],
                     key=lambda r: r['log_xi'])
    assert [r['log_xi'] for r in records] == [1, 2, 3]
    mpl.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.labelsize': 10,
        'axes.titlesize': 11, 'axes.linewidth': .7, 'xtick.labelsize': 8,
        'ytick.labelsize': 8, 'xtick.direction': 'in', 'ytick.direction': 'in',
        'xtick.top': True, 'ytick.right': True, 'legend.frameon': False,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'savefig.facecolor': 'white',
    })
    fig, axes = plt.subplots(1, 3, figsize=(8.2, 3.2))
    fig.subplots_adjust(left=.09, right=.986, bottom=.18, top=.80, wspace=.26)
    provenance = {
        'reference': PAPER, 'anchors_keV': anchors.tolist(),
        'upper_axis_scale': 'linear',
        'normalization': 'Fhat(E) = F_E(E) / ((F_E(E_low)+F_E(E_high))/2)',
        'layout': 'One row, three columns, with full-spectrum ratio insets.',
        'continuum_subtracted': False,
        'main_panels': 'Native energy samples, joined linearly; no extra Gaussian smoothing.',
        'ratio': 'Fhat_model / Fhat_DAO, after flux-preserving averaging onto common bins; full spectra including continuum.',
        'ratio_grid': 'Geometric-midpoint bin edges around native reflionx energy samples.',
        'ratio_limits_keV': [6.15, 7.10],
        'ratio_mask': 'Bins with non-finite or non-positive averaged DAO flux are excluded.',
        'runs': [],
    }
    for col, record in enumerate(records):
        xi = record['log_xi']
        models = load_models(record)
        prepared = {}
        info = {'log_xi': xi, 'DAO_hash': record['DAO_hash'],
                'DAO_iteration': record['DAO_iteration'], 'models': {}}
        top = axes[col]
        for name in NAMES:
            e, f = models[name]
            endpoint_flux = interp(e, f, anchors)
            assert np.all(np.isfinite(endpoint_flux)) and np.all(endpoint_flux > 0)
            scale = float(np.mean(endpoint_flux))
            a = endpoint_flux / scale
            assert np.isclose(a.mean(), 1, rtol=0, atol=2e-15)
            keep = (e > anchors[0]) & (e < anchors[1])
            plot_e = np.r_[anchors[0], e[keep], anchors[1]]
            normalized = interp(e, f, plot_e) / scale
            assert np.all(np.isfinite(normalized))
            top.plot(plot_e, normalized, color=COLORS[name], ls=STYLES[name], lw=1.15, label=name)
            top.plot(anchors, a, marker='o', ms=2.5, mfc='white', color=COLORS[name], ls='none')
            prepared[name] = {'e': e, 'f': f / scale, 'anchors': a}
            np.savetxt(HERE / 'data' / f'logxi{xi}_{name}.dat',
                       np.column_stack([plot_e, normalized]),
                       header='E_keV Fhat_E; complete spectrum, no continuum subtraction')
            info['models'][name] = {'endpoint_fluxes_input_units': endpoint_flux.tolist(),
                                    'normalization_divisor': scale,
                                    'normalized_anchor_mean': float(a.mean())}
        er = models['reflionx'][0]
        indices = np.where((er >= 6.15) & (er <= 7.10))[0]
        first, last = indices[0], indices[-1]
        edges = np.sqrt(er[first-1:last+1] * er[first:last+2])
        centers = er[indices]
        binned = {}
        for name in NAMES:
            model = prepared[name]
            averaged = bin_average(model['e'], model['f'], edges)
            check_e = np.r_[edges[0], model['e'][(model['e'] > edges[0]) & (model['e'] < edges[-1])], edges[-1]]
            before = TRAPEZOID(interp(model['e'], model['f'], check_e), check_e)
            after = np.sum(averaged * np.diff(edges))
            assert np.isclose(before, after, rtol=2e-13, atol=1e-14)
            info['models'][name]['rebin_integrated_flux_relative_error'] = float(abs(after-before)/before)
            binned[name] = averaged
        denominator = binned['DAO']
        valid = np.isfinite(denominator) & (denominator > 0)
        ratios = {name: np.divide(binned[name], denominator,
                                  out=np.full_like(denominator, np.nan), where=valid)
                  for name in NAMES[1:]}
        inset_position = [.62, .56, .34, .32] if xi < 3 else [.14, .56, .34, .32]
        inset = top.inset_axes(inset_position)
        inset.patch.set_facecolor('white')
        inset.patch.set_alpha(.97)
        inset.axhline(1, color=COLORS['DAO'], lw=.8)
        for name in NAMES[1:]:
            inset.plot(centers, ratios[name], color=COLORS[name], ls=STYLES[name], lw=.9,
                       marker='.', ms=2.7)
        inset.set_xlim(6.15, 7.10)
        inset.set_xticks([6.4, 6.7, 7.0])
        inset.yaxis.set_major_locator(MaxNLocator(3, min_n_ticks=2))
        inset.tick_params(labelsize=6.6, length=2, pad=1.5)
        inset.set_title('Model / DAO', fontsize=8, pad=2)
        inset.set_xlabel('Energy (keV)', fontsize=6.5, labelpad=1)
        np.savetxt(HERE / 'data' / f'logxi{xi}_ratios.dat',
                   np.column_stack([centers, edges[:-1], edges[1:], binned['DAO'],
                                    binned['reflionx'], binned['xillvercp'],
                                    ratios['reflionx'], ratios['xillvercp'], valid.astype(int)]),
                   header='E_keV bin_low_keV bin_high_keV Fhat_DAO Fhat_reflionx Fhat_xillvercp reflionx_over_DAO xillvercp_over_DAO valid\n'
                          'Ratios of complete normalized spectra; no continuum subtraction.')
        top.set_yscale('linear')
        top.yaxis.set_major_formatter(ScalarFormatter())
        top.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
        all_f = [interp(p['e'], p['f'], np.r_[anchors[0], p['e'][(p['e']>anchors[0]) & (p['e']<anchors[1])], anchors[1]])
                 for p in prepared.values()]
        ymax = max(float(v.max()) for v in all_f)
        top.set_ylim(0, ymax*1.15)
        top.yaxis.set_major_locator(MaxNLocator(5, min_n_ticks=3))
        top.set_title(rf'$\log\xi={xi}$', pad=9)
        top.set_xlabel('Energy (keV)')
        top.set_xlim(*anchors)
        top.set_xticks(np.arange(np.ceil(anchors[0]), anchors[1]+.1, 1))
        top.xaxis.set_minor_locator(AutoMinorLocator(5))
        top.tick_params(which='major', length=3.2, width=.65)
        top.tick_params(which='minor', length=1.8, width=.5)
        label_x = .04 if xi < 3 else .96
        top.text(label_x, .94, chr(ord('a')+col), transform=top.transAxes,
                 ha='left' if label_x < .5 else 'right', va='top', fontweight='bold', fontsize=9)
        files = [INPUTS / record[k] for k in ['DAO_file', 'DAO_params', 'reflionx_file', 'xillver_file']]
        info['source_sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
        expected = json.loads((INPUTS / 'plot_inputs.json').read_text())['runs']
        expected = next(r['input_sha256'] for r in expected if r['DAO_hash'] == record['DAO_hash'])
        assert info['source_sha256'] == expected, 'Benchmark source data changed'
        info['ratio_valid_bins'] = int(valid.sum())
        info['ratio_total_bins'] = len(valid)
        provenance['runs'].append(info)
        print(f'logxi={xi}: normalized anchors verified; {valid.sum()}/{len(valid)} ratio bins valid')
    axes[0].set_ylabel(r'$\widehat{F}_E$')
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(.54, .995),
               ncol=3, fontsize=10, handlelength=2.6, columnspacing=2.0)
    fig.savefig(output / 'pdf' / 'iron_line_zoom_logxi123.pdf')
    fig.savefig(output / 'iron_line_zoom_logxi123.png', dpi=600)
    (HERE / 'plot_manifest.json').write_text(json.dumps(provenance, indent=2)+'\n')
    if args.show:
        plt.show()
    plt.close(fig)


if __name__ == '__main__':
    main()
