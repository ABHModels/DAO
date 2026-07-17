# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) - [SemVer](https://semver.org/).

## [1.1.0] - 2026-07-05
- **Line-escape and thermal-feedback update for the DAO v1.0 model series.**

### Added
- `source/cloudy_interface.h`, `source/cloudy_interface_v2.cpp` - Add per-transition line records so bound-bound line power and opacity can be stored over the full slab before applying escape.
- `source/cloudy_interface_v2.cpp` - Add a second-pass, two-sided slab escape treatment for bound-bound lines, including line escape, electron-scattering escape, collisional quenching, and continuum destruction.
- `source/radiation.h`, `source/radiation.cpp`, `source/cloudy_interface_v2.cpp` - Add `rad.line_heat` to store continuum-destroyed line power and feed it back to Cloudy thermal balance through `hextra` on the following iteration.
- `source/save_results.cpp` - Add `line_heat_Pdest` and `heating_plus_new_line_heat` columns to `profile_iter*.dat` for thermal-feedback diagnostics.
- `source/cloudy_interface_v2.cpp`, `paper_data/figure_escape_probability/plot_escape_probability.py` - Add optional per-line escape diagnostics (`line_escape_lines_iter*.dat`, `line_escape_selected_iter*.dat`) and a plotting workflow based on converged model output.

### Changed
- Bound-bound line emissivity is now collected from Cloudy as optically thin local line power and deposited into the RT emissivity grid only after the full-column escape probability is evaluated.
- The implemented line survival factor is now `P = (beta + Pelec) * (1 + y) / (beta + Pelec + y + Pdest)`, where `beta` is two-sided line escape, `Pelec` is electron-scattering escape from the line core, `y = C_ul/A_ul` is collisional quenching, and `Pdest` is continuum destruction during trapping.
- H I Ly-alpha line power is converted back from Cloudy's internal escape treatment before DAO applies its slab-scale escape branching.
- `source/production.cpp` - Store line records for every depth cell during the Cloudy loop, call `apply_line_escape()` after the full column is available, and increase the outer-iteration cap from 50 to 100.

### Fixed
- Improve the physical treatment of bound-bound line escape. Previous versions attenuated line photons only with the redistribution-dependent escape probability `beta`. In v1.1.0, DAO uses the full escape-probability branching on the slab optical-depth grid: photons may escape through the line channel, be shifted out of the resonance core by electron scattering, be collisionally quenched, or be destroyed by continuum absorption. Only the continuum-destroyed line power is returned to Cloudy as an extra local heating term. Collisional quenching is not added as a separate DAO heating source because it is already part of Cloudy's local atomic and thermal balance; adding it again would double count that energy exchange.

## [1.0.0] - 2026-06-26
- **The first formal version of DAO**

### Changed
- `source/radiation.cpp:141` - Revert to our original definition $$\xi = (4\pi)^2 J / n_h$$, which follows [Tarter et al. 1969](https://ui.adsabs.harvard.edu/abs/1969ApJ...156..943T/abstract).

### Fixed
- `source/cloudy_interface_v2.cpp` (`extract_cloudy_output()`) - Subtract Cloudy's bound-electron Compton recoil opacity from `kabs`. Cloudy folds it into `opacity_abs` as absorption (the "No scattering opacity" command does not remove it), but our solver already handles all electron scattering via `ksct` and the kernel. Leaving it double-counted scattering as absorption, halving the 20-100 keV albedo (~0.95 -> ~0.5) and erasing the Compton hump.
- `source/cloudy_interface_v2.cpp` (`issue_depth()`, `issue_depth_lastest()`) - Pass the incident SED to Cloudy as the per-eV mean intensity `J0` (the F_nu shape), not `J0*wid/E` - the old factor was the log-grid spacing Delta(ln E) and distorted the shape by the grid resolution. Also floor SED values at `1e-30` instead of dropping them, preserving absorption troughs that Cloudy's log-log interpolation would otherwise bridge over.

## [0.1.3 beta] - 2026-06-18

### Changed

- `source/compton_cross_section.cpp:255` - Always use the total electron density to calculate the Compton-scattering cross section. The original version used only free electrons for Compton scattering, which caused temperature instability when the ionization parameter is low.
- `source/radiation.cpp:141` - Use the electron density to normalize the total flux.
  Original: 

  $$
  \xi = (4\pi)^2 J / n_h
  $$


  Now: 
  $$
  n_h \to n_e = 1.21 n_h
  $$
- `source/source.cpp` - Same as item 1, and removes some redundant code; also removes the 1e6 K limit on the electron temperature. At low ionization, this temperature approximation makes the iron Kα shoulder disappear. However, at low temperature (e.g. 1e4 K) the Compton redistribution function needs very high resolution to stay smooth; with the default resolution provided, jagged edges in the high-energy range are inevitable.

## [0.1.2 beta] - 2026-05-22

### Added

- `CloudyInput::issue_depth_lastest()` — post-convergence diagnostic run; dumps converged SED and iron ionization fractions.
- `source/dump_kernel_slice.cpp` — kernel slice dumper for debugging.
- `resolution/`, `paper_data/` — resolution-study scripts and paper figure data.

### Changed

- Outer-loop convergence (`production.cpp`): use **max** instead of mean of |Δlog T|, |Δlog ξ|; tolerance 1e-3 → 3e-3; max iters 15 → 20. Loop now stops only when **both** T and ξ are converged.
- `check_convergence()` args renamed `mean_dT/dXi` → `max_dT/dXi`.
- `source.cpp`: clamp kernel temperature to 1e6 K when `T_K < 1e6 K` (workaround for low-T numerical instability).

### Fixed

- `extract_cloudy_output()`: `esca0k2(tau)` is now applied only when redistribution is `ipCRD` (previously applied to all non-Lyα lines).

### Removed

- `save_cloudy_opacity()`, `save_line_labels()` — superseded by `save_results()`.
- Dead/commented blocks in `production.cpp` and `cloudy_interface_v2.cpp`.

## [0.1.0 beta] - 2026-05-08

- Initial Cloudy + Compton-RT coupling, redistribution kernel, short-characteristics solver.
