"""Reproduce both four-column native Cloudy/XSTAR comparisons from curated data.

No model execution, FITS reader, absolute paths, or compiled XSTAR required.
"""
from pathlib import Path
import argparse
import json
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/dao-cloudy-xstar-paper-mpl')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
BASE = Path(__file__).resolve().parent


def draw(density):
    case = BASE / f'lognH{density}'
    metadata = json.loads((case/'provenance.json').read_text())
    data = np.loadtxt(case/'data/native_comparison.dat')
    e = data[:, 0]*1000
    ct, xt, cr, xr = data[:, 1:5].T
    assert np.isfinite(data).all()
    assert np.allclose(data[:, 5], xt/ct)
    assert np.allclose(data[:, 6], xr/cr)
    anchors = [.55, .75]
    cr_norm = cr/np.interp(anchors, e/1000, cr).mean()
    xr_norm = xr/np.interp(anchors, e/1000, xr).mean()
    assert np.isclose(np.interp(anchors, e/1000, cr_norm).mean(), 1)
    assert np.isclose(np.interp(anchors, e/1000, xr_norm).mean(), 1)
    zoom = (e >= 550)&(e <= 750)
    saved_zoom = np.loadtxt(case/'data/oxygen_zoom.dat')
    calculated_zoom = np.c_[e[zoom]/1000, cr_norm[zoom], xr_norm[zoom],
                           xr_norm[zoom]/cr_norm[zoom], ct[zoom], xt[zoom], xt[zoom]/ct[zoom]]
    assert np.allclose(saved_zoom, calculated_zoom, rtol=1e-12, atol=1e-14)
    tc_front, tc_back = metadata['Cloudy_temperature_front_back_K']
    tx_front, tx_back = metadata['XSTAR_temperature_front_back_K']
    OUT = case/'figures'
    OUT.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size': 12, 'axes.linewidth': 1,
                         'xtick.direction': 'in', 'ytick.direction': 'in',
                         'svg.fonttype': 'none'})
    fig = plt.figure(figsize=(17, 5.1))
    grid = fig.add_gridspec(1, 4, wspace=.34)
    colors = ['#0072B2', '#D55E00']
    panels = [('Transmittance', ct, xt, (.2, 10)),
              ('Reflectance', cr, xr, (.2, 10)),
              ('Normalized reflectance', cr_norm, xr_norm, (.55, .75)),
              ('Transmittance', ct, xt, (.55, .75))]
    for j, (label, cy, xy, limits) in enumerate(panels):
        sub = grid[j].subgridspec(2, 1, height_ratios=[3, 1], hspace=.06)
        ax = fig.add_subplot(sub[0]); ratio_ax = fig.add_subplot(sub[1], sharex=ax)
        m = (e/1000>=limits[0])&(e/1000<=limits[1])
        for y, name, color, ls in [(cy, 'Cloudy', colors[0], '-'), (xy, 'XSTAR', colors[1], '--')]:
            ax.plot(e[m]/1000, y[m], color=color, ls=ls, lw=1.3, label=name)
        ratio_ax.plot(e[m]/1000, xy[m]/cy[m], color='#6e3c88', lw=.8)
        ratio_ax.axhline(1, color='0.45', ls=':', lw=.9)
        ax.set_ylabel(label)
        ratio_ax.set_ylabel('XSTAR/Cloudy', fontsize=10)
        ratio_ax.set_xlabel('Energy (keV)')
        ax.tick_params(labelbottom=False)
        for a in (ax, ratio_ax):
            a.set_xlim(*limits); a.tick_params(which='both', top=True, right=True)
        if j<2:
            ax.set_xscale('log')
        if j==0:
            ax.legend(frameon=False, loc='lower left')
            ax.text(.04, .96, rf'$n_{{\rm H}}=10^{{{density}}}\,\mathrm{{cm}}^{{-3}}$',
                    transform=ax.transAxes, va='top', fontsize=11)
            text = f'$T_{{\\rm back}}$\nCloudy: {tc_back/1e6:.3f} MK\nXSTAR: {tx_back/1e6:.3f} MK'
            ax.text(.04, .40, text, transform=ax.transAxes, va='top', fontsize=11)
            ax.set_ylim(0.6, 1.15)
        elif j in (1, 2):
            ax.set_yscale('log'); ratio_ax.set_yscale('log')
            if j==1:
                text = f'Cloudy: $T_{{\\rm front}}={tc_front/1e6:.3f}$ MK\nXSTAR: $T_{{\\rm front}}={tx_front/1e6:.3f}$ MK'
                ax.set_ylim(top=10)
                ax.text(.04, .96, text, transform=ax.transAxes, va='top', fontsize=11)
            else:
                ax.text(.96, .96, r'O VIII Ly$\alpha$', transform=ax.transAxes,
                        ha='right', va='top', fontsize=11)
        else:
            ax.text(.96, .96, r'O VIII Ly$\alpha$', transform=ax.transAxes,
                    ha='right', va='top', fontsize=11)
            low, high = min(cy[m].min(), xy[m].min()), max(cy[m].max(), xy[m].max())
            pad = max((high-low)*.15, 1e-3)
            ax.set_ylim(low-pad, high+pad)
    fig.subplots_adjust(left=.05, right=.99, bottom=.12, top=.98)
    fig.savefig(OUT/'native_comparison.png', dpi=230)
    fig.savefig(OUT/'native_comparison.svg')
    fig.savefig(OUT/'native_comparison.pdf')
    plt.close(fig)
    print(OUT/'native_comparison.pdf')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--density', type=int, choices=(15,18), help='Default: both densities')
    args = parser.parse_args()
    for density in ([args.density] if args.density else (15,18)):
        draw(density)
