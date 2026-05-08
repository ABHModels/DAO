# DAOv2.0 — X-ray Reflection Spectroscopy Model

**Disk-corona Atmosphere with Opacity**

A Compton scattering radiative transfer (RT) solver for computing X-ray reflection spectra from accretion disk atmospheres. The code couples **Cloudy** (Gary Ferland's spectral synthesis code) for atomic physics with custom second-order short characteristics RT, and uses **HEASoft/Xspec** models for corona spectra.

**Author:** Yimin Huang

---

## Table of Contents

- [Physics Overview](#physics-overview)
- [Architecture](#architecture)
- [Dependencies](#dependencies)
- [Environment Setup](#environment-setup)
- [Build](#build)
- [Usage](#usage)
- [Web UI](#web-ui)
- [Code Structure](#code-structure)
- [Algorithms & References](#algorithms--references)
- [Output Files](#output-files)
- [Important Notes](#important-notes)
- [License & Acknowledgements](#license--acknowledgements)

---

## Physics Overview

DAOv2.0 solves the angle- and energy-dependent radiative transfer equation in a plane-parallel slab (accretion disk atmosphere) illuminated by an X-ray corona from above and a thermal disk from below. The code iterates between:

1. **Cloudy** — computes emissivity j_nu, absorption opacity kappa_abs, and thermal structure at each depth point given the local radiation field
2. **RT solver** — solves the transfer equation with Compton scattering using second-order short characteristics and Lambda iteration

The Compton scattering kernel is computed exactly (not in the diffusion or Fokker-Planck approximation), following the relativistic formalism of Madej et al. (2017). The original Fortran source code for the exact redistribution function was generously shared by **Prof. Jerzy Madej** (University of Warsaw), whose pioneering work on exact Comptonisation made this project possible. We are deeply grateful for his generosity and his extraordinary contributions to the field. The Fortran code has been implemented in C++ for this project.

---

## Architecture

### Execution Flow

```
CLI arguments
    │
    ▼
read_params()  →  ModelParams (validate corona model + required params)
    │
    ▼
Bootstrap Cloudy  →  get energy grid (~3300 bins from Cloudy's rfield)
    │
    ▼
init_rt_grids()  →  angle (4-pt GL quadrature), depth (Thomson tau), energy grids
    │
    ▼
rad.illum.compute()  →  corona spectrum (via dispatch) + disk blackbody, normalise to F_x
    │
    ▼
Precompute KernelCache + ScatteringCache (load from disk or compute + save)
    │
    ▼
Dispatch:
    ├── test mode  →  run_test_rt()  (synthetic uniform slab, no Cloudy)
    └── production →  Outer loop:
                        ├── Cloudy depth sweep → extract j_nu, kappa_abs, kappa_sct
                        ├── RT solve (formal solution + Lambda iteration)
                        ├── Compute moments J, ionisation parameter xi
                        ├── Check convergence (mean |dT/T|, mean |d(log xi)|)
                        └── Save results every 5 iterations + iteration 1
```

### Corona Model Dispatch

| Model | Function | Parameters | Interface |
|-------|----------|------------|-----------|
| `powerlaw` | E^(1-Gamma) * exp(-E_lo/E) | -Gamma | Analytic |
| `cutoffpl` | E^(1-Gamma) * exp(-E/E_cut) * exp(-E_lo/E) | -Gamma, -Ecut | Analytic |
| `nthcomp` | Thermal Comptonisation | -Gamma, -kT_e, -kT_bb | Xspec `donthcomp_()` |
| `comptt` | Comptonisation | -kT_e, -kT_bb, -taup | Xspec `C_compTT()` |
| `blackbody` | Planck function | -kT_bb | Analytic |

### Disk-cached Precomputation

| Cache | Directory | File pattern | Contents |
|-------|-----------|-------------|----------|
| Compton kernel | `kernel/` | `kernel_norm_NE{}_NA{}_NT{}.bin` | Banded K(x,mu;x1,mu1,T), ~10% fill |
| Scattering sigma | `comp_cs/` | `scache_NE{}_NT{}.bin` | sigma(E,T) on 50-point T grid |

Both encode grid dimensions in the filename so different grids (test vs production) coexist. On first run they compute and save; subsequent runs load instantly.

---

## Dependencies

| Software | Version | Purpose |
|----------|---------|---------|
| **Cloudy** | C23+ | Atomic physics, emissivity/opacity, thermal equilibrium |
| **HEASoft** | 6.33+ | Xspec model libraries (nthcomp, comptt) |
| **g++** | C++17 | Compiler |
| **Python 3** | 3.8+ | Web UI (optional) |
| **Flask** | 2.0+ | Web UI server (optional, `pip install flask`) |

---

## Environment Setup

### 1. Cloudy

Build Cloudy as a library (`libcloudy.a`). The Makefile expects:
```
CLOUDY_SRC = /path/to/cloudy/source
CLOUDY_LIB = /path/to/cloudy/source
```
Edit these paths in `Makefile` to match your installation.

### 2. HEASoft / Xspec

**Required for nthcomp, comptt, and all production runs** (Xspec model libraries are linked at compile time).

```bash
export HEADAS=/path/to/heasoft-6.33.2/<arch>
source $HEADAS/headas-init.sh
```

**This must be sourced before every run** — the Xspec shared libraries need the `HEADAS` environment.

### 3. Python (optional, for Web UI)

```bash
pip install flask
```

---

## Build

```bash
# Edit Makefile paths for CLOUDY_SRC, CLOUDY_LIB, XSPEC_LIB if needed

make                    # build maindaocl (main executable)
make clean              # remove object files and binaries
make test_kernel        # standalone kernel test (no Cloudy dependency)
make test_kernel_fine   # fine-grained kernel validation (A23 sum rule)
make normalize_kernel   # normalize unnormalized kernel cache
```

---

## Usage

The `-corona` flag is **required**. Each model requires specific parameters:

```bash
# Cutoff power law
./maindaocl -corona cutoffpl -Gamma 2.0 -Ecut 300 -nh 16 -zeta 4 -frac 0.5

# Thermal Comptonisation (nthcomp)
./maindaocl -corona nthcomp -Gamma 2.0 -kT_e 60 -kT_bb 0.1 -nh 15 -zeta 3

# CompTT
./maindaocl -corona comptt -kT_e 50 -kT_bb 0.05 -taup 1.0 -Afe 3.0

# Test mode (no Cloudy, uniform slab)
./maindaocl -test_rt -corona nthcomp -Gamma 2.0 -kT_e 100 -kT_bb 0.05
```

### Full Parameter List

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `-corona` | string | **required** | Corona model: `powerlaw`, `cutoffpl`, `nthcomp`, `comptt`, `blackbody` |
| `-Gamma` | float | — | Photon index (powerlaw, cutoffpl, nthcomp) |
| `-Ecut` | float | — | High-energy cutoff [keV] (cutoffpl) |
| `-E_low_cut` | float | 0.1 | Low-energy exp cutoff [keV] (powerlaw, cutoffpl) |
| `-kT_e` | float | — | Electron temperature [keV] (nthcomp, comptt) |
| `-kT_bb` | float | — | Seed photon temperature [keV] (nthcomp, comptt, blackbody) |
| `-taup` | float | — | Plasma optical depth (comptt) |
| `-nh` | float | 15 | log hydrogen density [cm^-3] |
| `-zeta` | float | 3.0 | Ionisation parameter exponent: xi = 10^zeta |
| `-frac` | float | 100 | Flux ratio F_corona / F_disk |
| `-incidence` | float | 0.7071 | cos(theta) incidence angle (snapped to GL node) |
| `-Afe` | float | 1.0 | Iron abundance [solar] |
| `-kT_disk` | float | 0.35 | Disk blackbody temperature [eV] |
| `-test_rt` | flag | off | Enable test mode (skip Cloudy) |
| `-T_test` | float | 1e8 | Slab temperature [K] in test mode |

### Run Management

Results are saved to `results/<hash>/` where `<hash>` is an 8-character FNV-1a hash of all physics parameters. **Same parameters always produce the same hash** — re-running overwrites the same directory. Each run directory contains a `params.json` with the full parameter set.

---

## Web UI

A browser-based configurator and results viewer:

```bash
python3 ui.py
# Opens http://127.0.0.1:5200
```

### Pages

| Page | URL | Description |
|------|-----|-------------|
| **Configurator** | `/` | Parameter sliders, corona model selector, command generator, batch queue |
| **Results Viewer** | `/plots` | Select a run, view emergent spectra, mean intensity, temperature profile |
| **Parameter Reference** | `/docs` | Full documentation of all parameters with formulas |

The configurator generates commands — copy and run in your terminal. The results viewer reads `results/*/params.json` and plots data from the `.dat` files.

---

## Code Structure

```
DAOv2/
├── maindaocl.cpp              Main entry point
├── Makefile                   Build configuration
├── ui.py                      Web UI (Flask)
│
├── source/                    All C++ source modules
│   ├── compton_kernel.h/cpp           Exact azimuth-integrated kernel (Madej+ 2017)
│   ├── compton_cross_section.h/cpp    Klein-Nishina + relativistic sigma (Poutanen & Svensson 1996)
│   ├── compton_rt.h/cpp               Second-order short characteristics RT solver
│   ├── source.h/cpp                   Source function: thermal + Compton scattering
│   ├── radiation.h/cpp                RadField arrays, moments, convergence check
│   ├── corona_models.h/cpp            Corona spectrum dispatch (powerlaw, nthcomp, comptt, ...)
│   ├── cloudy_interface.h             Emissivity/opacity extraction interface
│   ├── cloudy_interface_v2.cpp        Cloudy emissivity/opacity extraction, normalization
│   ├── cloudy_exception.h             Cloudy error handling
│   ├── rt_grids.h/cpp                 Angle (GL quadrature), depth (Thomson tau), energy grids
│   ├── constants.h                    Physical constants (CGS), unit conversions
│   ├── params.h/cpp                   CLI parsing, parameter validation, hash computation
│   ├── save_results.h/cpp             Data output to results/<hash>/
│   ├── production.h/cpp               Cloudy-RT outer iteration loop
│   ├── test_rt.h/cpp                  Synthetic slab test mode
│   ├── normalize_kernel.cpp           Standalone tool: normalize cached kernel files
│   ├── test_kernel_norm.cpp           Standalone test: A23 normalization sum rule
│   ├── test_kernel_sym.cpp            Standalone test: kernel symmetry check
│   └── test_kernel_compare.cpp        Standalone test: compare against Madej Fortran code
│
├── plot/                      Plotting scripts (matplotlib)
├── image/                     UI assets
├── kernel/                    Cached Compton kernels (binary, generated at runtime)
├── comp_cs/                   Cached scattering cross-sections (binary, generated at runtime)
└── results/                   Output: results/<hash>/{emergent,moments,profile}_iter*.dat
```

---

## Algorithms & References

### Cloudy — Spectral Synthesis & Atomic Physics

- **Ferland G.J., Chatzikos M., Guzman F. et al., 2017, RMxAA, 53, 385** — "The 2017 Release of Cloudy"
  [ADS](https://ui.adsabs.harvard.edu/abs/2017RMxAA..53..385F) | [Cloudy website](https://gitlab.nublado.org/cloudy/cloudy)

Cloudy provides the atomic physics backend: emissivity, opacity, thermal equilibrium, and ionisation balance at each depth point. DAOv2.0 calls Cloudy as a subroutine via `cdInit()` / `cdDrive()`.

### Compton Scattering Kernel

The exact azimuth-integrated redistribution kernel K(x, mu; x1, mu1, T):

- **Madej J., Rozanska A., Majczyna A., Nalezyta M., 2017, MNRAS, 469, 2032** — Exact redistribution function `profil()`, differential cross section
  [ADS](https://ui.adsabs.harvard.edu/abs/2017MNRAS.469.2032M)
- **Nagirner D.I., Poutanen J., 1993, A&A, 275, 325** — Azimuth-integrated kernel formalism
  [ADS](https://ui.adsabs.harvard.edu/abs/1993A%26A...275..325N)
- **Poutanen J., Svensson R., 1996, ApJ, 470, 249** — Normalization condition (Eq. A23), exact cross section
  [ADS](https://ui.adsabs.harvard.edu/abs/1996ApJ...470..249P)

Implementation: 32-point Gauss-Laguerre quadrature over electron Lorentz factor, 6-point Gauss-Legendre for azimuthal integration.

### Compton Cross Section

Relativistic scattering cross section sigma(E, T) averaged over a Maxwellian electron distribution:

- **Poutanen J., Svensson R., 1996, ApJ, 470, 249** — Exact relativistic cross section
  [ADS](https://ui.adsabs.harvard.edu/abs/1996ApJ...470..249P)
- Three regimes: Thomson (low T, low E), Klein-Nishina (low T, high E), exact relativistic (high T)
- Modified Bessel function K2 via **Abramowitz & Stegun, 1972, Handbook of Mathematical Functions** polynomial approximation

### Radiative Transfer Solver

Second-order short characteristics with Lambda iteration:

- **Hubeny I., Mihalas D., 2015, Theory of Stellar Atmospheres, Princeton University Press, Section 12.4** — Second-order short characteristics formal solution
  [Publisher](https://press.princeton.edu/books/hardcover/9780691163291/theory-of-stellar-atmospheres)
- **Suleimanov V., Poutanen J., Werner K., 2012, A&A, 545, A120** — Application to neutron star spectrum
  [ADS](https://ui.adsabs.harvard.edu/abs/2012A%26A...545A.120S)

### Emissivity Extraction from Cloudy

- **Garcia J., Kallman T.R., 2010, ApJ, 718, 695** — Emissivity normalization (Eq. 12), X-ray reflection modelling
  [ADS](https://ui.adsabs.harvard.edu/abs/2010ApJ...718..695G)
- Conversion: `j_nu = (ConEmitLocal[1] + DiffuseLineEmission) * anu(j) * eV_to_erg / widflx(j)`
- Opacity: `kappa_abs = opacity_abs + OpacStatic`, `kappa_sct = n_e * sigma_compton(E,T)`

### Corona Models (HEASoft/Xspec)

- **Zdziarski A.A., Johnson W.N., Magdziarz P., 1996, MNRAS, 283, 193** — nthcomp thermal Comptonisation
  [ADS](https://ui.adsabs.harvard.edu/abs/1996MNRAS.283..193Z)
- **Titarchuk L., 1994, ApJ, 434, 570** — comptt Comptonisation model
  [ADS](https://ui.adsabs.harvard.edu/abs/1994ApJ...434..570T)
- **Arnaud K.A., 1996, ASP Conf. Ser., 101, 17** — XSPEC spectral fitting package
  [ADS](https://ui.adsabs.harvard.edu/abs/1996ASPC..101...17A) | [HEASoft](https://heasarc.gsfc.nasa.gov/docs/software/heasoft/)

### Ionisation Parameter

- **Tarter C.B., Tucker W.H., Salpeter E.E., 1969, ApJ, 156, 943** — Definition of ionisation parameter xi = 4pi F_x / n_H
  [ADS](https://ui.adsabs.harvard.edu/abs/1969ApJ...156..943T)

---

## Output Files

Each run saves to `results/<hash>/` with:

| File | Contents |
|------|----------|
| `params.json` | Full parameter set (JSON) |
| `emergent_iter{NNN}.dat` | Emergent intensity at surface: E, I_corona, I_disk, I(mu) for each angle |
| `moments_iter{NNN}.dat` | Angular moments J0, J2, J3 at all depths and energies |
| `profile_iter{NNN}.dat` | Depth profile: tau, T, n_e, heating, cooling, log_xi |

### Data Units

| Array | Units |
|-------|-------|
| Energy E | eV |
| j_nu (emissivity) | erg cm^-3 s^-1 eV^-1 |
| kappa_abs, kappa_sct | cm^-1 |
| I_nu (specific intensity) | erg cm^-2 s^-1 eV^-1 sr^-1 |
| I_corona, I_disk | erg cm^-2 s^-1 eV^-1 |

---

## Important Notes

1. **Cloudy headers**: Files that include Cloudy headers must include `cddefines.h` before any other Cloudy header.

2. **fopen**: Cloudy redefines `fopen` to a compile error. Use `open_data()` for Cloudy I/O, or `#undef fopen` before raw file operations (e.g., before `cdInit()`).

3. **SED files**: `table SED` input files must floor values at `1e-30` — Cloudy's log-interpolation asserts on zeros.

4. **Energy range**: The Compton kernel is initialized once covering 0.1–1000 keV. Below 0.1 keV, Thomson scattering applies.

5. **Temperature cache**: Constants `T_CACHE_LO` (10^4 K), `T_CACHE_HI` (10^9 K), `N_T_CACHE` (50) in `constants.h` are shared by `KernelCache` and `ScatteringCache`.

6. **ConEmitLocal index**: Use `ConEmitLocal[1]` (iteration phase), not `[0]` (search phase, always zero).

7. **Unit convention**: All radiation quantities use **per-eV** units internally. Use `phys::eV_to_erg` (not `EN1RYD`) for conversions.

8. **Run deduplication**: Same parameters always produce the same hash and overwrite the same output directory.

---

## License & Acknowledgements

This project relies on and gratefully acknowledges:

- **Cloudy** by Gary Ferland and collaborators — atomic physics and spectral synthesis backend.
- **HEASoft / Xspec** by NASA HEASARC — Comptonisation model libraries (`nthcomp`, `comptt`).
- **Prof. Jerzy Madej** (University of Warsaw) — for sharing the original Fortran source code of the exact Compton redistribution function, on which this implementation is based.

If you use DAOv2.0 in published work, please cite the references listed in the [Algorithms & References](#algorithms--references) section, in particular Madej et al. (2017) for the Compton kernel and Ferland et al. (2017) for Cloudy.

For questions, suggestions, or collaboration, please open an issue on GitHub or contact the author.
