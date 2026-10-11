# Figure 3 — DAO vs `reflionx` vs `xillvercp`

**Author:** Yimin Huang · Fudan University; University of Bristol · huangym23@m.fudan.edu.cn

Updated on 2026-10-09 using the latest converged DAO spectra. The original
three-panel layout, model colours, reference spectra, and 20–50 keV
normalisation are retained. The panels compare spectral *shape*.

Shared input parameters: `nthcomp`, Γ = 2, kT_e = 60 keV,
kT_bb = 0.01 keV, log n_H = 15, and solar iron abundance.
The DAO runs have `frac=-1` (no incident disk component), 48 depth cells,
angle-averaged scattering, and the coupled thermal/transfer solver.

The directory retains its historical `figure_benchmark_xillvernorm` name.
The replacement DAO runs use **`reflionx_norm=1`**, as recorded in their
parameter files. This boundary-flux convention is separate from the
20–50 keV scaling applied to every plotted curve.

## Input snapshots

| log ξ | DAO hash | Final iteration | Reference-file hash |
|---|---|---|---|
| 1 | `9e14af6b` | 63 | `1e4178a5` |
| 2 | `effac594` | 35 | `603b2ef4` |
| 3 | `a518ab75` | 33 | `feb16067` |

- `dao_<DAO hash>.dat`: snapshot of the final upper-face specific intensities.
- `pa<DAO hash>.json`: parameters for that DAO snapshot.
- `spectra_<reference hash>.dat`: unchanged `reflionx` reference spectrum.
- `xillver_<reference hash>.dat`: unchanged `xillvercp` reference spectrum.
- `benchmark_inputs.json`: source hashes, iteration numbers, convergence metadata,
  and the reference parameters used to check consistency with DAO. These
  parameters are retained from the original reference snapshots; the deleted
  old DAO parameter files are not required for plotting.
- `plot_inputs.json`: plotted inputs, SHA-256 hashes, and normalisation factors.
- `compare_reflionx_xi_scan.{png,pdf}`: updated publication figures.

Only the current DAO snapshots are retained. The REFLIONX and XILLVER-CP
reference files keep their original identifiers and remain active inputs.

## Flux calculation

The DAO reflected flux is the outgoing normal flux at the illuminated face,
`F_E = 2*pi*sum(w_i*mu_i*I_i)` over `mu_i>0`. Exact eight-point
Gauss–Legendre weights are used; header nodes are checked for consistency.
The old plotting function used the outgoing contribution to mean intensity;
the updated function computes the normal flux used by the current benchmark.

The requested DAO incidence cosine is 0.7071. On its eight-direction grid,
the actual incoming direction has `|mu_inc| = 0.7966664774`, which is the
value annotated in the figure. Each curve is then divided by its own mean
F_E over 20–50 keV, after interpolating the references onto the DAO grid.

## Run
```bash
python plot_compare_reflionx.py
```
All required input files are local to this folder. The command saves the
PDF and 600 dpi PNG without opening a window. For interactive zoom and pan:

```bash
python plot_compare_reflionx.py --show
```
