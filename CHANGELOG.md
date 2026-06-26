# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) · [SemVer](https://semver.org/).

## [1.0.0] - 2026-06-26
- **The first formal version of DAO**

### Changed
- `source/radiation.cpp:141` - Revert to our original definition $$\xi = (4\pi)^2 J / n_h$$, which follows [Tarter et al. 1969](https://ui.adsabs.harvard.edu/abs/1969ApJ...156..943T/abstract).

### Fixed
- `source/cloudy_interface_v2.cpp` (`extract_cloudy_output()`) - Subtract Cloudy's bound-electron Compton recoil opacity from `kabs`. Cloudy folds it into `opacity_abs` as absorption (the "No scattering opacity" command does not remove it), but our solver already handles all electron scattering via `ksct` and the kernel. Leaving it double-counted scattering as absorption, halving the 20–100 keV albedo (~0.95 → ~0.5) and erasing the Compton hump. 
- `source/cloudy_interface_v2.cpp` (`issue_depth()`, `issue_depth_lastest()`) - Pass the incident SED to Cloudy as the per-eV mean intensity `J0` (the F_ν shape), not `J0·wid/E` — the old factor was the log-grid spacing Δ(ln E) and distorted the shape by the grid resolution. Also floor SED values at `1e-30` instead of dropping them, preserving absorption troughs that Cloudy's log-log interpolation would otherwise bridge over.

## [0.1.3 beta] - 2026-06-18

### Changed
- `source/compton_cross_section.cpp:255` - Always use the total electron density to calculate the Compton-scattering cross section. The original version used only free electrons for Compton scattering, which caused temperature instability when the ionization parameter is low.
- `source/radiation.cpp:141` - Use the electron density to normalize the total flux. 
    Original: $$\xi = (4\pi)^2 J / n_h$$ 
    Now: $$n_h \to n_e = 1.21 n_h$$
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
