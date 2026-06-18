# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) · [SemVer](https://semver.org/).

## [2.2.1] - 2026-06-18
- `source/compton_cross_section.cpp:255` - Always use total electron densities to calculate the Compton scattering cross section. The original one only use the free electron for compton scattering, it'll cause the temperature instability when ionization parameter is low.
- `source/radiation.cpp:141` - Use electron density to normalize the total flux. 
    Original: $$ \xi = \frac{(4\pi)^2J/n\_h}$$ 
    Now : $$n\_h \to n\_e = 1.21n\_h$$
- `source/source.cpp` - same as item 1 and remove some redundant code; remove 1e6 K limit for electron temperature, when ionization is low, this approximation temperature will make the iron Kalpha shoulder disappear. However, at low temperature (e.g. 1e4K) the compton redistribution function need very high resoltion to keep smooth. If one use default resolution we provided, there will inevitably have jagged edges in the high-energy range.


## [2.2.0] - 2026-05-22

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

## [2.1.0] - 2026-05-08

- Initial Cloudy + Compton-RT coupling, redistribution kernel, short-characteristics solver.
