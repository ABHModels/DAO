#!/usr/bin/env python
"""Schematic of the corona--slab geometry (publication style).

Single-column figure for the paper: slab band along x-axis, corona above
indicated by a subtle ellipse, incident ray (solid + dashed continuation
through the slab) and emergent ray, with angles theta_inc and theta_view
measured from the surface normal (the y-axis).
"""

import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Circle, Ellipse
from matplotlib.colors import LinearSegmentedColormap

# ── Nature/ApJ-style rcParams (small, sans-serif, clean math) ───────────────
mpl.rcParams.update({
    'font.family':         'sans-serif',
    'font.sans-serif':     ['Helvetica', 'Arial', 'DejaVu Sans'],
    'mathtext.fontset':    'dejavusans',
    'font.size':            8,
    'axes.linewidth':       0.5,
    'lines.linewidth':      0.9,
    'pdf.fonttype':         42,
    'ps.fonttype':          42,
})

# ── Geometry (illustrative angles, change at will) ──────────────────────────
THETA_INC_DEG = 55          # incidence angle from y-axis (surface normal)
THETA_EMI_DEG = 36          # emergent (viewing) angle from y-axis
RAY_LEN       = 2.5         # solid arrow length
EXT_LEN       = 1.35        # dashed extension length (through the slab)

# Subdued palette
C_INC  = '#1f6fb5'          # incident   - muted blue
C_EMI  = '#b32a2a'          # emergent   - muted red
C_AX   = '0.25'             # axes and neutral text - dark grey
C_COR  = '#9ec3e0'          # corona ellipse - pale blue
C_SLAB = '#a8a8a8'          # slab hatch line colour


# ── helper: wavy "photon" line from p0 to p1, with an arrowhead at p1 ───────
def wavy_photon(ax, p0, p1, color, n_waves=6, amp=0.07, lw=1.1, zorder=4):
    p0 = np.asarray(p0, dtype=float)
    p1 = np.asarray(p1, dtype=float)
    d = p1 - p0
    L = np.linalg.norm(d)
    if L == 0:
        return
    u    = d / L
    perp = np.array([-u[1], u[0]])

    # Wavy body: stop slightly short to leave room for the arrowhead.
    t_end = 0.93
    t = np.linspace(0.0, t_end, 260)
    base = p0[None, :] + np.outer(t, d)
    # sin-envelope tapers amplitude to 0 at both ends -> clean attachment.
    envelope = np.sin(np.pi * t / t_end)
    disp = (amp * envelope * np.sin(2 * np.pi * n_waves * t / t_end))[:, None] \
           * perp[None, :]
    xy = base + disp
    ax.plot(xy[:, 0], xy[:, 1], color=color, lw=lw, zorder=zorder,
            solid_capstyle='round', solid_joinstyle='round')

    # Arrowhead from end-of-wave (no displacement there) to the true endpoint.
    arrow_start = base[-1]
    ax.annotate('', xy=p1, xytext=arrow_start,
                arrowprops=dict(arrowstyle='-|>', color=color, lw=lw,
                                mutation_scale=9),
                zorder=zorder)


def main():
    fig, ax = plt.subplots(figsize=(3.6, 3.2))

    theta_inc = np.deg2rad(THETA_INC_DEG)
    theta_emi = np.deg2rad(THETA_EMI_DEG)
    d_inc = np.array([np.sin(theta_inc), -np.cos(theta_inc)])   # down-right
    d_emi = np.array([np.sin(theta_emi),  np.cos(theta_emi)])   # up-right

    # ── slab band (subtle hatching below x-axis) ────────────────────────────
    slab_xlim = [-3.2, 3.2]
    ax.fill_between(slab_xlim, -1.05, 0.0, facecolor='none',
                    edgecolor=C_SLAB, hatch='///', lw=0.0, zorder=0)
    # Thin solid line along the slab surface (the x-axis itself)
    ax.plot(slab_xlim, [0, 0], color=C_AX, lw=0.7, zorder=2)
    ax.text(3.05, -0.52, r'slab', ha='right', va='center',
            fontsize=8.5, color=C_AX, style='italic', zorder=3)

    # ── corona: smooth elliptical Gaussian blue haze + bold label ──────────
    # Built as an imshow of a 2D *anisotropic* Gaussian, using a custom
    # colormap that fades from fully transparent at the edges to a deep blue
    # at the centre — same look as slides/allgeometry.png (flat horizontal
    # ellipse, wider than tall, not a square blob).
    # Shift slightly to +x so the left edge of the imshow box doesn't get
    # clipped by the axis at x = -3.20: with this corona_cx the box left edge
    # sits exactly at x = -3.20, where the Gaussian is already below the 0.25
    # transparency cutoff, so it fades naturally instead of being chopped.
    corona_cx, corona_cy = -1.55, 2.00
    w_box, h_box = 3.3, 1.30           # plot-coordinate extent of the haze

    nx, ny = 420, 200
    xx = np.linspace(-1.0, 1.0, nx)
    yy = np.linspace(-1.0, 1.0, ny)
    Xc, Yc = np.meshgrid(xx, yy)
    sigma_x, sigma_y = 0.60, 0.32      # anisotropic -> ellipse (aspect ~2:1)
    gauss = np.exp(-(Xc ** 2 / (2.0 * sigma_x ** 2)
                     + Yc ** 2 / (2.0 * sigma_y ** 2)))

    # Colormap stays *fully transparent* below value 0.25 so the box edges
    # (where the Gaussian is still ~0.2) fade to invisible, removing the
    # spurious rectangular halo.
    cmap_corona = LinearSegmentedColormap.from_list(
        'corona',
        [(0.00, (0.55, 0.72, 0.88, 0.00)),   # transparent
         (0.25, (0.55, 0.72, 0.88, 0.00)),   # still transparent
         (0.50, (0.50, 0.68, 0.86, 0.25)),   # soft haze
         (0.80, (0.35, 0.58, 0.82, 0.62)),   # mid blue
         (1.00, (0.20, 0.43, 0.70, 0.88))],  # deep core
        N=256,
    )

    extent = [corona_cx - w_box / 2, corona_cx + w_box / 2,
              corona_cy - h_box / 2, corona_cy + h_box / 2]
    ax.imshow(gauss, extent=extent, origin='lower', cmap=cmap_corona,
              interpolation='bicubic', aspect='auto',
              vmin=0.0, vmax=1.0, zorder=1)

    # Bold uppercase label, with a subtitle, in the convention of the
    # reference figure (slides/allgeometry.png).
    ax.text(corona_cx, corona_cy + 0.08, 'CORONA',
            ha='center', va='center', fontsize=10.5, fontweight='bold',
            color='#163152', zorder=5)
    ax.text(corona_cx, corona_cy - 0.13, '(X-ray source)',
            ha='center', va='center', fontsize=7, fontstyle='italic',
            color='#1a3a55', zorder=5)

    # ── x and y axes (thin grey arrows; y-axis IS the surface normal) ───────
    arrow_ax = dict(arrowstyle='-|>', color=C_AX, lw=0.55,
                    mutation_scale=7)
    ax.annotate('', xy=(3.15, 0), xytext=(-3.2, 0), arrowprops=arrow_ax)
    ax.annotate('', xy=(0, 3.05), xytext=(0, -1.35), arrowprops=arrow_ax)
    ax.text(3.15, -0.22, r'$x$', fontsize=9, color=C_AX, style='italic')
    ax.text(-0.20, 3.00, r'$y\equiv\hat n$', fontsize=9, color=C_AX,
            ha='right', va='top')

    # ── incident photon (wavy line in, dashed straight continuation in slab)
    upstream = -d_inc * RAY_LEN
    wavy_photon(ax, upstream, (0.0, 0.0), color=C_INC, n_waves=6, amp=0.075,
                lw=1.05, zorder=4)
    ext_end = d_inc * EXT_LEN
    ax.plot([0, ext_end[0]], [0, ext_end[1]], color=C_INC, lw=0.8,
            ls=(0, (4, 2)), zorder=4)

    # ── emergent photon (wavy line out) ─────────────────────────────────────
    emi_end = d_emi * RAY_LEN
    wavy_photon(ax, (0.0, 0.0), emi_end, color=C_EMI, n_waves=6, amp=0.075,
                lw=1.05, zorder=4)

    # ── angle arcs (dashed) and labels (italic Greek) ───────────────────────
    # Arcs are drawn dashed and slightly extended toward the wavy ray, but
    # their radius stays well below RAY_LEN so they never reach the photon.
    arc_ls = (0, (3, 2))                  # custom dash pattern for the arcs
    arc_r_inc = 1.00                      # slight extension toward incident ray
    ax.add_patch(Arc((0, 0), 2 * arc_r_inc, 2 * arc_r_inc,
                     theta1=90, theta2=90 + THETA_INC_DEG,
                     color=C_INC, lw=0.7, ls=arc_ls, zorder=3))
    a_inc = np.deg2rad(90 + THETA_INC_DEG / 2)
    r_lab_inc = 1.30
    ax.text(r_lab_inc * np.cos(a_inc), r_lab_inc * np.sin(a_inc),
            r'$\theta_{\rm inc}$', ha='center', va='center',
            color=C_INC, fontsize=8.5)

    arc_r_emi = 0.95                      # slight extension toward emergent ray
    ax.add_patch(Arc((0, 0), 2 * arc_r_emi, 2 * arc_r_emi,
                     theta1=90 - THETA_EMI_DEG, theta2=90,
                     color=C_EMI, lw=0.7, ls=arc_ls, zorder=3))
    a_emi = np.deg2rad(90 - THETA_EMI_DEG / 2)
    r_lab_emi = 1.25
    ax.text(r_lab_emi * np.cos(a_emi), r_lab_emi * np.sin(a_emi),
            r'$\theta_{\rm view}$', ha='center', va='center',
            color=C_EMI, fontsize=8.5)

    # ── ray labels (small italic, perpendicular offset OUTSIDE each ray) ────
    # Push each label to the side of its ray that is away from the y-axis,
    # so the angle labels theta_inc / theta_view sit cleanly on the y-axis
    # side and the ray labels don't collide with them.
    perp_inc = np.array([ d_inc[1], -d_inc[0]])    # outward (away from +y)
    perp_emi = np.array([ d_emi[1], -d_emi[0]])    # outward (away from +y)
    inc_mid = upstream * 0.55 + perp_inc * 0.20
    emi_mid = emi_end  * 0.55 + perp_emi * 0.22
    ax.text(inc_mid[0], inc_mid[1], 'incident', ha='center', va='center',
            color=C_INC, fontsize=8, style='italic',
            rotation=-(90 - THETA_INC_DEG), rotation_mode='anchor')
    ax.text(emi_mid[0], emi_mid[1], 'emergent', ha='center', va='center',
            color=C_EMI, fontsize=8, style='italic',
            rotation=(90 - THETA_EMI_DEG), rotation_mode='anchor')

    # ── cleanup: no spines, no ticks ─────────────────────────────────────────
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ['top', 'right', 'bottom', 'left']:
        ax.spines[s].set_visible(False)

    ax.set_xlim(-3.2, 3.2)
    ax.set_ylim(-1.35, 3.10)
    ax.set_aspect('equal')

    plt.tight_layout(pad=0.1)

    plt.show()
    out_dir = '/Users/ym.huang/cloudy_test/image'
    os.makedirs(out_dir, exist_ok=True)
    out_png = f'{out_dir}/incidence_geometry.png'
    out_pdf = f'{out_dir}/incidence_geometry.pdf'
    fig.savefig(out_png, dpi=600, bbox_inches='tight')
    fig.savefig(out_pdf, bbox_inches='tight')
    print(f'Saved: {out_png}')
    print(f'Saved: {out_pdf}')


if __name__ == '__main__':
    main()
