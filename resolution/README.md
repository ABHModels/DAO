# Resolution Files for Compton RT Code

Energy-dependent resolving-power profiles (`R(E) = E/ΔE`) used to define
the spectral grid of the Compton RT solver. The generator script builds
a smooth Gaussian "hump" in `log10(E)` so that resolution peaks around
the Compton hump / Fe Kα region and falls back to a coarser baseline at
low and high energies.

## Author

- **Name:** Yimin Huang
- **Affiliation:** Fudan University; University of Bristol
- **Email:** huangym23@m.fudan.edu.cn

## Contents

| File                  | Description |
|-----------------------|-------------|
| `make_smooth_hump.py` | Generator for `smooth_hump.ini`; also plots `R(E)` and `1/R(E)`. |
| `smooth_hump.ini`     | Resolving-power table consumed by the RT code (Ryd vs. R). |
| `make_log_resolution.py` | Generates `log_resolution.ini` with constant logarithmic spacing; no plot. |
| `log_resolution.ini` | Constant-resolution alternative for Cloudy's continuum mesh. |

## Simple logarithmic mesh

From the DAO repository root:

```bash
python3 resolution/make_log_resolution.py --resolution 300
```

This writes only `resolution/log_resolution.ini`. Defaults are 1 eV to
1000 keV; change them with `--emin-ev` and `--emax-kev`. Use `--output`
to choose a different file. The default output is next to the script,
regardless of the working directory.

Cloudy uses a constant nominal `delta(ln E) = 1/R` in this interval,
adjusting the bin count to fit its boundaries. Thus larger `R` gives finer
spacing. Cloudy may also insert special atomic edges. Outside this interval,
the file retains the existing low/high-energy settings (10 and 33.333333).
DAO's ready-to-use mesh is supplied as `config/continuum_mesh.ini`.
To generate a custom mesh there, use `--output config/continuum_mesh.ini`.
Set `CLOUDY_DATA_PATH` as described in the main [README](../README.md);
Cloudy's installed data files do not need to be replaced.

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

For questions, please contact Yimin Huang (Fudan University; University of Bristol) —
huangym23@m.fudan.edu.cn.
