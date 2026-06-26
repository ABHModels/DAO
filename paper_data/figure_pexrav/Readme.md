# Figure 4 — Incidence-angle sensitivity of the reflection spectrum

**Author:** Yimin Huang · Fudan University · huangym23@m.fudan.edu.cn (hyimin0924@gmail.com)

Two panels, both compared against `pexrav` (reflection-only):
- **(a)** DAO incidence-averaged emergent spectrum (at μ_obs = 0.7).
- **(b)** DAO reflected spectrum at each illumination-angle bin.

All runs share: `cutoffpl` corona, Γ = 2, E_cut = 300 keV, log n_H = 15,
log ξ = 0, kT_disk = 0.35 keV.

## Files
- `emergent_<hash>.dat` — final-iteration emergent spectrum (E, incident corona,
  incident disk, then I at 8 GL μ nodes). `38ea5924_inc_avg` is the
  incidence-averaged run (panel a); the eight other hashes are per-incidence
  runs (panel b).
- `params_<hash>.json` — run parameters (the `incidence` field drives panel b).
- `pexrav.dat` — `pexrav` reflection-only reference (E [keV], ΔE/2, value).
- `plot_compare_pexrav_inc.py`, `compare_pexrav_inc.{png,pdf}`.

The solver snaps each requested `-incidence` to the nearest Gauss–Legendre
node, so the eight per-incidence runs collapse onto 4 distinct angle bins
(0.1834, 0.5255, 0.7967, 0.9603); the script does this grouping automatically.

## Run
```bash
python plot_compare_pexrav_inc.py
```
No arguments needed; all paths are local to this folder.
