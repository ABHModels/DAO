# Figure 1 — Escape probability for subordinate lines

**Author:** Yimin Huang · Fudan University; University of Bristol · huangym23@m.fudan.edu.cn

Left panel: schematic of the cell-by-cell solver. Right panel: escape
probability β_ℓ vs line-center optical depth τ_ℓ. `beta_K2` and `beta_PRD`
are exact ports of Cloudy's `esca0k2` and `esc_PRD_1side` (`rt_escprob.cpp`);
the constant-`p_w` dashed curves are a schematic CRD+wing simplification
(Hummer 1982; Ferland et al. 2017), not a port of `esc_CRDwing_1side`.

## Files
- `plot_escape_probability.py` — self-contained; computes everything analytically.
- `escape_probability.pdf` — the figure.

## Run
```bash
python plot_escape_probability.py
```
No data files or arguments needed.
