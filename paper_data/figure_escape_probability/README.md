# Figure — Line escape in DAO: escape-probability branching

**Author:** Yimin Huang · Fudan University; University of Bristol · huangym23@m.fudan.edu.cn

Two-row small-multiple figure. Four columns = four representative lines chosen
by DAO's `write_line_escape_diagnostics()`, one per escape-physics regime:
(a) O VIII Lyα 18.97 Å (strong escape), (b) O VIII 16.01 Å (continuum
destruction), (c) Ne VIII 770.4 Å (strong trapping), (d) O VII 21.81 Å (weak
attenuation).

- **Top row** — WHAT it does: local line power vs Thomson depth τ_T, comparing
  the optically-thin emissivity with the power that actually escapes; the shaded
  band is the power removed by trapping/destruction.
- **Bottom row** — HOW it is done: the survival factor
  `P = (β + P_el)(1 + y) / (β + P_el + y + P_dest)` decomposed into its four
  competing channels — line (Sobolev) escape β, electron-scattering escape
  P_el, continuum destruction P_dest, and collisional quenching y = C_ul/A_ul.

Data are from DAO model `603b2ef4`, latest iteration (025).

## Files
- `line_escape_lines_603b2ef4.dat` — per-line slab summary (rank, ip, energy,
  slab-integrated thin/escaped power, mean channel ratios, label).
- `line_escape_selected_603b2ef4.dat` — per-depth channel profiles for the four
  selected representative lines.
- `plot_escape_probability.py` — reads the two `.dat` files above; all paths are
  local to this folder.
- `escape_probability.{pdf,png}` — the figure.
- `escape_probability_caption.tex` — LaTeX caption.

## Run
```bash
python plot_escape_probability.py
```
No arguments needed; all paths are local to this folder.
