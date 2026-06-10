# Figure 8 — Angle-averaged vs directional Compton kernel

**Author:** Yimin Huang · Fudan University · huangym23@m.fudan.edu.cn (hyimin0924@gmail.com)

Two DAOv2.0 cutoffpl runs with identical physics (Γ = 2, E_cut = 300 keV,
n_H = 10^15 cm^-3, log ξ = 2, incidence cos θ = 0.70 → snapped GL node
μ_inc = 0.797). Only the Compton kernel differs:

- `a1f9c7b5` — directional kernel (`angsca = true`)
- `1e89bf78` — angle-averaged kernel (`angsca = false`)

2×2 small multiples, one cell per viewing GL node
(μ = 0.183, 0.526, 0.797, 0.960). Each cell stacks a spectrum panel
(raw emergent I_E, solid = directional, dashed = angle-averaged, gap shaded)
over a ratio panel (angle-averaged / directional, dotted y = 1 parity line).

## Files
- `angavg_vs_directional.py` — reads the two `.dat` files + `params.json` and writes the figure.
- `emergent_a1f9c7b5_iter020.dat` — emergent I(μ,E), directional kernel (final iteration).
- `emergent_1e89bf78_iter020.dat` — emergent I(μ,E), angle-averaged kernel (final iteration).
- `params.json` — run parameters (from `a1f9c7b5`; `1e89bf78` is identical except `angsca`).
- `angavg_vs_directional.{png,pdf}` — the figure.

## Run
```bash
python angavg_vs_directional.py
```
No arguments needed.
