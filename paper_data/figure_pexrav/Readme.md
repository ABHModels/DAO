# Figure 4 — DAO / PEXRAV reflection comparison

**Author:** Yimin Huang · Fudan University; University of Bristol · huangym23@m.fudan.edu.cn

The two-panel figure compares the final DAO emergent spectra at viewing cosine
μ_view = 0.7 with the reflection-only reference in `pexrav.dat`.

- **(a)** Single-angle illumination at μ_inc = 0.796666 with the
  angle-averaged Compton kernel (dashed orange).
- **(b)** Isotropic illumination and four single-angle incident beams, all
  with the angle-dependent Compton kernel.

Every DAO run uses a cutoff power law, Γ = 2, E_cut = 300 keV,
log(n_H/cm⁻³) = 15, log ξ = 0, kT_disk = 0.35 keV, A_Fe = 1,
and the same energy grid. The requested beam cosines 0.18, 0.52, 0.79, and
0.96 snap to the positive Gauss–Legendre nodes 0.183435, 0.525532,
0.796666, and 0.960290.
The script reads the outgoing μ nodes from each spectrum and interpolates
linearly between nodes to obtain μ_view = 0.7.

| Run hash | Final spectrum in `results/` | Incidence | Compton kernel |
| --- | --- | --- | --- |
| `11168045` | `emergent_iter014.dat` | isotropic (`incidence=-2`) | angle dependent |
| `310ae399` | `emergent_iter013.dat` | requested μ=0.18 | angle dependent |
| `04b092df` | `emergent_iter012.dat` | requested μ=0.52 | angle dependent |
| `3b9a421c` | `emergent_iter011.dat` | requested μ=0.79 | angle dependent |
| `607c5c2f` | `emergent_iter012.dat` | requested μ=0.96 | angle dependent |
| `3c9a43af` | `emergent_iter011.dat` | requested μ=0.79 | angle averaged |

The `emergent_<hash>.dat` and `params_<hash>.json` files here are copies of
those final spectra and run parameters, so the plot is reproducible without
the runtime `results/` directories. The source spectrum files contain
energy in eV, incident corona and disk intensities, and emergent intensities
at eight angular nodes. `pexrav.dat` contains energy in keV, half-bin width,
and reflection-only intensity. PEXRAV is interpolated in log energy and log
intensity over its positive range; the figure does not extrapolate outside it.
Each plotted curve is separately divided by its mean intensity over
0.1–1000 keV.

Run `python plot_compare_pexrav_inc.py` in this directory to regenerate
`compare_pexrav_inc.png` and `compare_pexrav_inc.pdf`.
