# Figure 3 (xillver normalisation) — DAO vs `reflionx` vs `xillvercp`

**Author:** Yimin Huang · Fudan University; University of Bristol · huangym23@m.fudan.edu.cn

Emergent reflection spectra from DAO compared against `reflionx` and
`xillvercp` at log ξ = 1, 2, 3 (one panel each). All runs share: `nthcomp`
corona, Γ = 2, kT_e = 60 keV, kT_bb = 0.01 keV, kT_disk = 0.35 keV,
log n_H = 15, μ_inc = cos 45°. Each model is normalised to its mean over the
20–50 keV continuum, so the panels compare spectral *shape*.

DAO here uses the **`xillver` flux-normalisation convention**. The
companion folder `figure_benchmark/` is the same comparison with DAO's
native flux normalisation; the plotting script is identical and only the DAO
input spectra differ.

## Files (one set per run hash: `1e4178a5`=ξ1, `603b2ef4`=ξ2, `feb16067`=ξ3)
- `dao_<hash>.dat` — DAO emergent intensity (E, incident, I at 8 μ nodes).
- `spectra_<hash>.dat` — `reflionx` spectrum (E [keV], E·F_E).
- `xillver_<hash>.dat` — `xillvercp` spectrum.
- `pa<hash>.json` — run parameters.
- `plot_compare_reflionx.py`, `compare_reflionx_xi_scan.{png,pdf}`.

## Run
```bash
python plot_compare_reflionx.py
```
No arguments needed; all paths are local to this folder.
