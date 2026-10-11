# Figure 8 — Angle-averaged vs directional Compton kernel

**Author:** Yimin Huang · Fudan University; University of Bristol · huangym23@m.fudan.edu.cn

Two DAO cutoffpl runs with identical physics (Γ = 2, E_cut = 300 keV,
n_H = 10^15 cm^-3, log ξ = 2, incidence cos θ = 0.70 → snapped GL node
μ_inc = 0.797). Only the Compton kernel differs:

- `60d8ae0b` — directional kernel (`angsca = true`), final iteration 20
- `5fd8ac78` — angle-averaged kernel (`angsca = false`), final iteration 21

2×2 small multiples, one cell per viewing GL node
(μ = 0.183, 0.526, 0.797, 0.960). Each cell stacks a spectrum panel
(raw emergent I_E, solid = directional, dashed = angle-averaged, gap shaded)
over a ratio panel (angle-averaged / directional, dotted y = 1 parity line).

## Files
- `angavg_vs_directional.py` — reads the two `.dat` files + `params.json` and writes the figure.
- `emergent_60d8ae0b_iter020.dat` — emergent I(μ,E), directional kernel (final iteration).
- `emergent_5fd8ac78_iter021.dat` — emergent I(μ,E), angle-averaged kernel (final iteration).
- `params.json` and `params_5fd8ac78.json` — new run parameters.
- `angavg_vs_directional.{png,pdf}` — the figure.

Only the two current spectra listed above are retained.

## Run
```bash
python angavg_vs_directional.py
```
No arguments needed.
