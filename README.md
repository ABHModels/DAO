# DAO — X-ray Reflection Spectroscopy Model

**Disk-corona Atmosphere with Opacity**

A Compton scattering radiative transfer (RT) solver for computing X-ray reflection spectra from accretion disk atmospheres. The code couples **Cloudy** (Gary Ferland's spectral synthesis code) for atomic physics with a custom RT solver based on a cubic Bezier interpolant, and uses **HEASoft/Xspec** models for corona spectra.

**Author:** Yimin Huang (Fudan University; University of Bristol) · huangym23@m.fudan.edu.cn

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

DAO solves the angle- and energy-dependent radiative transfer equation in a plane-parallel slab (accretion disk atmosphere) illuminated by an X-ray corona from above and a thermal disk from below. The code iterates between:

1. **Cloudy** — computes emissivity j_nu, absorption opacity kappa_abs, and thermal structure at each depth point given the local radiation field
2. **RT solver** — solves the transfer equation with Compton scattering using a cubic Bezier interpolant for the formal solution, with Lambda iteration

The Compton scattering kernel is computed exactly (not in the diffusion or Fokker-Planck approximation), following the relativistic formalism of [Madej et al. (2017)](https://ui.adsabs.harvard.edu/abs/2017MNRAS.469.2032M). The original Fortran source code for the exact redistribution function was generously shared by **Prof. Jerzy Madej** (University of Warsaw), whose pioneering work on exact Comptonisation made this project possible. We are deeply grateful for his generosity and his extraordinary contributions to the field. The Fortran code has been implemented in C++ for this project.

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
rad.illum.compute()  →  corona spectrum (via dispatch) + disk blackbody, normalise to ionization parameter
    │
    ▼
Precompute KernelCache (load from disk or compute + save)
    │
    ▼
Dispatch:
    ├── test mode  →  run_test_rt()  (no Cloudy):
    │                   compps : isothermal pure-scattering slab illuminated from
    │                            bottom by a blackbody seed, benchmarked against
    │                            Xspec compPS (Poutanen & Svensson 1996)
    └── production →  Outer loop:
                        ├── Cloudy depth sweep → extract j_nu, kappa_abs, kappa_sct
                        ├── RT solve (cubic Bezier formal solution + Lambda iteration)
                        ├── Compute moments J, ionisation parameter xi
                        ├── Check convergence (max |dT/T|, max |d(log xi)|)
                        └── Save results
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
| Compton kernel | `$COMPTON_CACHE_DIR/kernel/` | `kernel_norm_NE{}_NI{}_NT{}.bin` | Banded K(x,mu;x1,mu1,T), ~10% fill |

The cache directory is controlled by the `COMPTON_CACHE_DIR` environment variable (defaults to the current working directory). Grid dimensions are encoded in the filename so different grids (test vs production) coexist. On first run the cache is computed and saved; subsequent runs load instantly.

---

## Dependencies

| Software | Version | Purpose |
|----------|---------|---------|
| **Cloudy** | C25+ (we also test for C23)| Atomic physics, emissivity/opacity, thermal equilibrium | 
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

**Custom continuum-mesh resolution file (required).** DAO ships with a tailored Cloudy resolution file at `resolution/smooth_hump.ini` that boosts the spectral resolving power around the iron K region (~30 keV, R ~ 600) while keeping it modest elsewhere. You must install it as Cloudy's `continuum_mesh.ini`:

```bash
# Back up Cloudy's default mesh first
cp /path/to/cloudy/data/continuum_mesh.ini /path/to/cloudy/data/continuum_mesh.ini.bak

# (Optional) Choose your own resolution before installing. Edit the parameters at
# the top of resolution/make_smooth_hump.py — R_base, R_peak, E_center, sigma,
# N_zones, and the energy bounds — then regenerate the table (the default value has good accuracy):
python resolution/make_smooth_hump.py        # writes resolution/smooth_hump.ini (+ a preview PNG)

# Install the DAOv2 resolution file
cp resolution/smooth_hump.ini /path/to/cloudy/data/continuum_mesh.ini
```

The resolving power follows a Gaussian "hump" in log energy,
`R(E) = R_base + (R_peak − R_base)·exp(−½ z²)` with `z = (log₁₀E_keV − log₁₀E_center)/sigma`,
peaked at `E_center` and falling to `R_base` elsewhere. See `resolution/README.md` for the parameters and file format.

**Without this step, the energy grid produced by Cloudy will be too coarse/large to run the scattering kernel correctly.**

### 2. HEASoft / Xspec

**Required for nthcomp, comptt, and all production runs** (Xspec model libraries are linked at compile time).

```bash
export HEADAS=/path/to/heasoft-6.33.2/<arch>
source $HEADAS/headas-init.sh
```

**This must be sourced before every run** — the Xspec shared libraries need the `HEADAS` environment.

### 3. Compton kernel cache directory

The Compton kernel is precomputed once and cached on disk under a `kernel/` subdirectory. By default the code looks for this in the current working directory; to use a shared location (recommended, since the kernel files are large and parameter-independent), set:

```bash
export COMPTON_CACHE_DIR=/path/to/shared/cache
# Kernel files will be looked up / written under $COMPTON_CACHE_DIR/kernel/
```

If `COMPTON_CACHE_DIR` is unset the code falls back to `./kernel/`. On first run the cache is computed and saved (this can take a while); every subsequent run loads it instantly.

### 4. Python (optional, for Web UI)

```bash
pip install flask
```

---

## Build

```bash
# Edit Makefile paths for CLOUDY_SRC, CLOUDY_LIB, XSPEC_LIB if needed

make                    # build maindaocl (main executable)
make clean              # remove object files and binaries
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

# Test mode — compPS benchmark slab (Poutanen & Svensson 1996): isothermal, pure
# scattering, illuminated by a bottom blackbody seed. kT_e is the slab temperature.
./maindaocl -test_rt compps -corona blackbody -kT_e 60 -kT_bb 0.1 -tau 0.5
```

### Full Parameter List

**Corona model (required)**

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `-corona` | string | **required** | Corona model: `powerlaw`, `cutoffpl`, `nthcomp`, `comptt`, `blackbody` |
| `-Gamma` | float | unset | Photon index Γ (powerlaw, cutoffpl, nthcomp) |
| `-Ecut` | float | unset | High-energy cutoff [keV] (cutoffpl) |
| `-E_low_cut` | float | 0.1 | Low-energy exponential cutoff [keV] (powerlaw, cutoffpl) |
| `-kT_e` | float | unset | Electron temperature [keV] (nthcomp, comptt) |
| `-kT_bb` | float | unset | Seed photon temperature [keV] (nthcomp, comptt, blackbody) |
| `-taup` | float | unset | Plasma optical depth (comptt) |

**Slab / illumination**

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `-nh` | float | 15 | log₁₀ hydrogen density [cm⁻³] |
| `-zeta` | float | 3.0 | Ionisation parameter exponent: ξ = 10^zeta |
| `-frac` | float | -1 | Flux ratio F_corona / F_disk; `≤ 0` (default) → corona only (no disk component) |
| `-incidence` | float | 0.7071 | cos(θ) of the corona incidence angle (snapped to the nearest GL node) |
| `-kT_disk` | float | 0.35 | Disk blackbody temperature [keV] (thermal component illuminating from below) |
| `-Afe` | float | 1.0 | Iron abundance [solar] |

**Solver / kernel**

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `-angsca` | bool | true | Scattering kernel: `true`/`1`/`yes` → angle-dependent `KernelCache`; `false`/`0`/`no` → angle-averaged `avgKernelCache` |

**Test mode (skip Cloudy)**

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `-test_rt <mode>` | flag + string | off | Enable test mode (no Cloudy). Mode token is **required**: `compps` (compPS benchmark slab) or `test_avg`. |
| `-kT_e` | float | unset | Slab uniform temperature [keV] in `compps` test mode (re-uses the corona `-kT_e` flag). |
| `-tau` | float | 0.5 | Slab vertical Thomson optical depth (`compps` test mode) |

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
│   ├── compton_rt.h/cpp               RT solver (cubic Bezier interpolant + Lambda iteration)
│   ├── source.h/cpp                   Source function: thermal + Compton scattering
│   ├── radiation.h/cpp                RadField arrays, moments, convergence check
│   ├── corona_models.h/cpp            Corona spectrum dispatch (powerlaw, nthcomp, comptt, ...)
│   ├── cloudy_interface.h             Emissivity/opacity extraction interface
│   ├── cloudy_interface_v2.cpp        Cloudy emissivity/opacity extraction
│   ├── cloudy_exception.h             Cloudy error handling
│   ├── rt_grids.h/cpp                 Angle (GL quadrature), depth (Thomson tau), energy grids
│   ├── constants.h                    Physical constants (CGS), unit conversions
│   ├── params.h/cpp                   CLI parsing, parameter validation, hash computation
│   ├── save_results.h/cpp             Data output to results/<hash>/
│   ├── production.h/cpp               Cloudy-RT outer iteration loop
│   └── test_rt.h/cpp                  Synthetic slab test mode
│
├── resolution/                Custom continuum-mesh files for Cloudy
│   └── smooth_hump.ini        Boosted resolving power around Fe K (install as Cloudy's continuum_mesh.ini)
├── plot/                      Plotting scripts (matplotlib)
├── image/                     UI assets
├── kernel/                    Cached Compton kernels (binary, generated at runtime)
└── results/                   Output: results/<hash>/{emergent,moments,profile}_iter*.dat
```

---

## Algorithms & References

### Cloudy — Spectral Synthesis & Atomic Physics

- **Gunasekera C.M., van Hoof P.A.M., Dehghanian M. et al., 2025, arXiv:2508.01102** — most recent release of Cloudy (C25)
  [ADS](https://ui.adsabs.harvard.edu/abs/2025arXiv250801102G) | [Cloudy website](https://gitlab.nublado.org/cloudy/cloudy)

Cloudy provides the atomic physics backend: emissivity, opacity, thermal equilibrium, and ionisation balance at each depth point. DAO calls Cloudy as a subroutine via `cdInit()` / `cdDrive()`.

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

Short characteristics formal solution using a cubic Bezier interpolant, with Lambda iteration:

- **Auer L., 2003, ASP Conf. Ser., 288, 3** — "Insertion of Cubic Bezier Splines into Short-Characteristics Solutions of the Radiative Transfer Equation"
  [ADS](https://ui.adsabs.harvard.edu/abs/2003ASPC..288....3A)
- **de la Cruz Rodriguez J., Piskunov N., 2013, ApJ, 764, 33** — Bezier formal solutions in radiative transfer
  [ADS](https://ui.adsabs.harvard.edu/abs/2013ApJ...764...33D)
- **Hubeny I., Mihalas D., 2015, Theory of Stellar Atmospheres, Princeton University Press, Section 12.4** — Short characteristics formal solution and Lambda iteration
  [Publisher](https://press.princeton.edu/books/hardcover/9780691163291/theory-of-stellar-atmospheres)
- **Suleimanov V., Poutanen J., Werner K., 2012, A&A, 545, A120** — Application to neutron star spectrum
  [ADS](https://ui.adsabs.harvard.edu/abs/2012A%26A...545A.120S)

### Corona Models (HEASoft/Xspec)

- **Zdziarski A.A., Johnson W.N., Magdziarz P., 1996, MNRAS, 283, 193** — nthcomp thermal Comptonisation
  [ADS](https://ui.adsabs.harvard.edu/abs/1996MNRAS.283..193Z)
- **Titarchuk L., 1994, ApJ, 434, 570** — comptt Comptonisation model
  [ADS](https://ui.adsabs.harvard.edu/abs/1994ApJ...434..570T)
- **Arnaud K.A., 1996, ASP Conf. Ser., 101, 17** — XSPEC spectral fitting package
  [ADS](https://ui.adsabs.harvard.edu/abs/1996ASPC..101...17A) | [HEASoft](https://heasarc.gsfc.nasa.gov/docs/software/heasoft/)

### Ionisation Parameter

- **Tarter C.B., Tucker W.H., Salpeter E.E., 1969, ApJ, 156, 943** — Definition of ionisation parameter xi = 4pi F_x / n_H, i.e., in the code, we calculate xi by (4pi)^2 J/n_H
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

## License & Acknowledgements

### License

The original DAO source code in this repository is released under the **MIT License** (see [`LICENSE`](LICENSE)).

DAO **depends on** the following third-party software, each of which carries its own separate license and must be obtained independently (DAO does not bundle their distributions):

| Component | License | Source |
|-----------|---------|--------|
| **Cloudy** (Gary J. Ferland and collaborators) — atomic physics / spectral synthesis backend | [zlib license](https://opensource.org/licenses/Zlib) | <https://gitlab.nublado.org/cloudy/cloudy> |
| **HEASoft / XSPEC** (NASA HEASARC) — Comptonisation model libraries (`nthcomp`, `comptt`) | NASA open-source (HEASARC) | <https://heasarc.gsfc.nasa.gov/docs/software/heasoft/> |
| **Compton redistribution function** — `source/compton_kernel.cpp` is a C++ port of the publicly available Fortran code by Jerzy Madej (please cite [Madej, Różańska, Majczyna & Należyta, 2017, MNRAS, 469, 2032](https://ui.adsabs.harvard.edu/abs/2017MNRAS.469.2032M)) | as published by the author | <https://www.astrouw.edu.pl/~jm/software.html> |

Cloudy's zlib license permits free use, modification, and redistribution. DAO uses Cloudy as an **unmodified** backend: it links against a separately installed Cloudy and accesses it only through Cloudy's public API — headers, library functions, and data structures. **DAO does not modify, copy, or redistribute Cloudy's source code**, so no altered Cloudy source is shipped. A few physics steps are handled in DAO's own code rather than delegated to Cloudy:

- the per-line escape probabilities — computed by **calling** Cloudy's `rt_escprob` routines;
- the removal of Cloudy's bound-electron Compton-recoil opacity term — **reproducing**, on DAO's own grid, the calculation in Cloudy's `opacity_addtotal.cpp`, because Cloudy folds that term into its opacity and exposes no API to subtract it;
- the inner-shell fluorescence emission — **reproducing** the calculation in Cloudy's `prt_lines.cpp`, using Cloudy's `t_yield` atomic data and inner-shell photoionization rates.

These reimplementations reference Cloudy's data structures and methods (cited in-line), but a physical formula plus the required API symbol names are not themselves Cloudy source code. Cloudy remains under its [zlib license](https://opensource.org/licenses/Zlib) (© 1978–2025 Gary J. Ferland and others).

### Acknowledgements

We gratefully acknowledge **Prof. Jerzy Madej** (University of Warsaw) for his publicly available Fortran implementation of the exact Compton redistribution function, on which our C++ port is based, and for his extraordinary contributions to the theory of exact Comptonisation.

If you use DAO in published work, please cite the references listed in the [Algorithms & References](#algorithms--references) section, in particular Madej et al. (2017) for the Compton kernel and Gunasekera et al. (2025) for Cloudy.

For questions, suggestions, or collaboration, please open an issue on GitHub or contact the author.
