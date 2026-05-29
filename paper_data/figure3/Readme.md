# Figure 3 — DAOv2.0 vs `reflionx` vs `xillvercp` over a log ξ scan

**Author:** Yimin Huang · Fudan University · huangym23@m.fudan.edu.cn (hyimin0924@gmail.com)

Emergent reflection spectra from DAOv2.0 compared against `reflionx` and
`xillvercp` at log ξ = 1, 2, 3 (one panel each). All runs share: `nthcomp`
corona, Γ = 2, kT_e = 60 keV, kT_bb = 0.01 keV, kT_disk = 0.35 keV,
log n_H = 15, μ_inc = cos 45°. Each model is normalised to its mean over the
20–50 keV continuum, so the panels compare spectral *shape*.

## Files (one set per run hash: `c63d0c48`=ξ1, `8b877447`=ξ2, `50bd6dd6`=ξ3)
- `dao_<hash>.dat` — DAOv2.0 emergent intensity (E, incident, I at 8 μ nodes).
- `spectra_<hash>.dat` — `reflionx` spectrum (E [keV], E·F_E).
- `xillver_<hash>.dat` — `xillvercp` spectrum.
- `pa<hash>.json` — run parameters.
- `plot_compare_reflionx.py`, `compare_reflionx_xi_scan.{png,pdf}`.

## Run
```bash
python plot_compare_reflionx.py
```
No arguments needed; all paths are local to this folder.
