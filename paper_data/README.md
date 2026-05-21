# Paper Data — Compton Scattering Radiative Transfer with Cloudy

This repository contains the data, source code, and plotting scripts
required to reproduce the figures in the paper.

## Author

- **Name:** Yimin Huang
- **Affiliation:** Fudan University
- **Email:** huangym23@m.fudan.edu.cn (alt: hyimin0924@gmail.com)

## Contents

Each subfolder corresponds to one figure in the paper and contains its
own README describing the files, how to (re)generate the data, and how
to reproduce the figure.

| Folder    | Figure | Topic |
|-----------|--------|-------|
| `figure1` | Fig. 1 | Escape probability for subordinate lines (CRD + Voigt wings); validated against Cloudy's `rt_escprob.cpp`. |
| `figure2` | Fig. 2 | 2-D slices of the Compton redistribution kernel `K(μ_out, E_out)` at fixed `(E_in, μ_in, T)` for Fe Kα and the Compton hump. |
| `figure3` | Fig. 3 | Emergent reflection spectra from DAOv2.0 compared with `reflionx` and `xillvercp` over a log ξ scan. |

## How to use

1. Enter the relevant subfolder.
2. Read its `README.md` / `Readme.md` for the figure-specific build and
   plotting instructions.
3. For C++ tools (e.g. `figure2/dump_kernel_slice.cpp`), the file must
   be copied into the main RT solver's `source/` directory and built
   against the Compton kernel routines — see the local README.
4. Python plotting scripts only depend on `numpy`, `matplotlib`, and the
   `.dat`/`.json` files distributed alongside them.

## License & citation

If you use the code or data here, please cite the corresponding paper
(see the paper bibliography for the exact reference) and contact the
author for any questions.

## Contact

For questions, please contact:

- Yimin Huang (Fudan University) — huangym23@m.fudan.edu.cn
