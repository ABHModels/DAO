# Paper Data

Data and self-contained plotting scripts that reproduce the figures in the DAO paper.

**Author:** Yimin Huang · Fudan University; University of Bristol · huangym23@m.fudan.edu.cn

## Figures

| Folder | Fig. | Topic |
|--------|------|-------|
| `figure_escape_probability` | 1 | Escape probability for subordinate lines (vs Cloudy `rt_escprob.cpp`). |
| `figure_angdep_CR` | 2 | 2-D slices of the Compton redistribution kernel (Fe Kα and the Compton hump). |
| `figure_benchmark` | 3 | DAO vs `reflionx` and `xillvercp` over a log ξ scan (native DAO normalisation). |
| `figure_benchmark_xillvernorm` | 3 | Same comparison, `xillver` flux-normalisation convention. |
| `figure_pexrav` | 4 | Incidence-angle sensitivity of the reflection spectrum vs `pexrav`. |
| `figure_compps` | 6 | Comptonised slab: DAO emergent spectra and limb darkening vs `compPS`. |
| `figure_compare_angdep` | 8 | Angle-averaged vs directional Compton kernel. |
| `figure_geometry` | — | Model geometry schematic. |

## Usage

Each folder is self-contained: enter it and run its script with no arguments
(e.g. `python plot_compare_reflionx.py`). See the folder's `Readme.md` for details.

- Most scripts need only `numpy` + `matplotlib`.
- `figure_compps` also needs **HEASoft / PyXspec** (it builds the `compPS` model).

## Citation

If you use this data or code, please cite the DAO paper.
