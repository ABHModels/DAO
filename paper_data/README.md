# Paper Data — Compton Scattering Radiative Transfer with Cloudy

Data and plotting scripts to reproduce the paper's figures.

**Author:** Yimin Huang · Fudan University · huangym23@m.fudan.edu.cn (alt: hyimin0924@gmail.com)

## Figures

| Folder    | Topic |
|-----------|-------|
| `figure1` | Escape probability for subordinate lines; validated against Cloudy `rt_escprob.cpp`. |
| `figure2` | 2-D slices of the Compton redistribution kernel for Fe Kα and the Compton hump. |
| `figure3` | DAOv2.0 reflection spectra vs `reflionx` and `xillvercp` over a log ξ scan. |
| `figure4` | Incidence-angle sensitivity of the reflection spectrum vs `pexrav`. |
| `figure5` | Iron ionization fraction vs Thomson optical depth (DAOv1.0 vs DAOv2.0). |
| `figure6` | Comptonized slab: DAOv2.0 emergent spectra and limb-darkening law vs `compPS`. |

## How to use

Each folder is self-contained: its data, plotting script, and figure all live
together, with one script per figure. Just enter a folder and run its script
with no arguments, e.g. `python plot_kernel_redist.py`. See each folder's
`Readme.md` for details.

- Most scripts need only `numpy` + `matplotlib`.
- `figure5` also needs `astropy` and `pandas`.
- `figure6` needs **HEASoft / PyXspec** (it builds the `compPS` model) — see its Readme.

## Citation

If you use this code or data, please cite the corresponding paper and contact
the author with any questions.
