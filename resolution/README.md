# Resolution Files for Compton RT Code

Energy-dependent resolving-power profiles (`R(E) = E/ΔE`) used to define
the spectral grid of the Compton RT solver. The generator script builds
a smooth Gaussian "hump" in `log10(E)` so that resolution peaks around
the Compton hump / Fe Kα region and falls back to a coarser baseline at
low and high energies.

## Author

- **Name:** Yimin Huang
- **Affiliation:** Fudan University
- **Email:** huangym23@m.fudan.edu.cn

## Contents

| File                  | Description |
|-----------------------|-------------|
| `make_smooth_hump.py` | Generator for `smooth_hump.ini`; also plots `R(E)` and `1/R(E)`. |
| `smooth_hump.ini`     | Resolving-power table consumed by the RT code (Ryd vs. R). |

## Resolving-power formula

```
R(E) = R_base + (R_peak - R_base) * exp(-0.5 * z^2)
z    = (log10(E_keV) - log10(E_center)) / sigma_dex
```

Edit the parameters at the top of `make_smooth_hump.py`
(`R_base`, `R_peak`, `E_center`, `sigma`, `N_zones`, energy bounds) and
re-run to regenerate the `.ini` table.

## Usage

```bash
python make_smooth_hump.py
```

Produces `smooth_hump.ini` and `smooth_hump_resolution.png`.

## File format (`smooth_hump.ini`)

- Line 1: magic number `10 08 08`.
- Comment lines start with `#`.
- Data lines: `upper_limit_Ryd   resolving_power`.
- Final line: `0   <closing_R>` marks the upper bound.

## Contact

For questions, please contact Yimin Huang (Fudan University) —
huangym23@m.fudan.edu.cn.
