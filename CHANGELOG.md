# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) · [SemVer](https://semver.org/).

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
