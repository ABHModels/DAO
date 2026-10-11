# Paper Data

Data and self-contained plotting scripts that reproduce the figures in the DAO paper.

**Author:** Yimin Huang · Fudan University; University of Bristol · huangym23@m.fudan.edu.cn

## Contents

| Directory | Contents |
|-----------|----------|
| `figure_benchmark_xillvernorm` | DAO / REFLIONX / XILLVER-CP at log ξ = 1, 2, 3; current DAO runs `9e14af6b`, `effac594`, `a518ab75`; final convergence values and LaTeX table. The historical directory name is retained, but the DAO inputs use REFLIONX-compatible incident-flux normalization. |
| `figure_benchmark_xillvernorm/iron_line_zoom` | Fe K zoom at the three ionizations, normalized at 5 and 8 keV, with reference/DAO ratio insets; no continuum subtraction or EW annotations. |
| `figure_all_survival_off_lognH18` | Broad-band and O VIII comparison at log nH = 18, log ξ = 3, Γ = 1.8; standard DAO versus all bound–bound survival factors set to one with the atmosphere re-equilibrated, alongside REFLIONX and XILLVER-CP. Includes experiment provenance. |
| `figure_pexrav` | Cold-reflection benchmark against PEXRAV, including isotropic and directional illumination. |
| `figure_compps` | Isothermal pure-scattering benchmark against compPS: spectra and angular dependence using the shared production transfer solver. |
| `figure_compare_angdep` | Angle-dependent versus angle-averaged Compton scattering at four viewing directions. |
| `figure_angdep_CR` | Compton redistribution-kernel slices for 6.4 and 40 keV incident photons. |
| `figure_escape_probability` | Local line emission, survival factors, and direct-escape/electron-scattering contributions for four He II transitions; run `64e9e118`, iteration 11. |
| `figure_fe_kalpha_production_logxi12` | Local Fe Kα production at log ξ = 1 and 2: all Fe, Fe I–V, and Fe XVII–XXII; includes extraction code and depth-integrated fractions. |
| `figure_fe_survival_logxi3` | Fixed-atmosphere Fe XXV w+x control at log ξ = 3: DAO, DAO with these two survival factors set to one, and XILLVER-CP. No EW annotations. |
| `figure_cloudy_xstar_nh15_nh18` | Native Cloudy/XSTAR at log nH = 15 and 18, log ξ = 3, NH = 10²¹ cm⁻²: transmittance, reflectance, O VIII zooms and ratios; curated plot data, complete model inputs and self-contained plotting code. Modified XSTAR executable/source excluded. |
| `figure_geometry` | Model-geometry schematic script. |

## Usage

Consult each directory's README for plotting commands and physical definitions.
Most saved-spectrum plots require only NumPy and Matplotlib. The Fe K zoom
reads the input snapshots in its parent benchmark directory.

- For compPS, use `python comppsGenerator.py --reuse-reference` to plot saved
  inputs; regenerating its reference spectra requires HEASoft / PyXspec.
- The escape-probability script also copies its output into the repository's
  `Figure/` directory.
- The geometry script saves its PNG and PDF beside the script, regardless of the working directory.
- Rebuilding the convergence table requires the original `results/` logs; the
  final table and extracted values are already included.
- Experiment drivers and emissivity-extraction tools may require the original
  model runs and Cloudy/HEASoft installation. Plotting saved inputs does not
  rerun those calculations.

## Citation

If you use this data or code, please cite the DAO paper.
