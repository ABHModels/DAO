# DAO — X-ray Reflection Spectroscopy Model

**Disk-corona Atmosphere with Opacity**

DAO computes rest-frame X-ray reflection spectra from accretion-disc atmospheres in X-ray binaries and AGN. Its distinguishing feature is an **exact, angle- and energy-dependent Compton-redistribution kernel** — fully relativistic quantum-electrodynamic scattering, not the angle-averaged Gaussian or Fokker–Planck approximations used by other reflection codes. DAO couples **Cloudy** (C25) for atomic physics, a **custom Compton radiative-transfer solver** (cubic Bézier short-characteristics formal solution with Lambda iteration), and **HEASoft/Xspec** corona spectra.

**Author:** Yimin Huang (Fudan University; University of Bristol) · huangym23@m.fudan.edu.cn
**Version:** 1.0.0 · MIT License

---

## Table of Contents

- [Quick Start](#quick-start)
- [Dependencies & Environment Setup](#dependencies--environment-setup)
- [Usage](#usage)
- [Outputs](#outputs)
- [Web UI](#web-ui)
- [How It Works](#how-it-works)
- [References, License & Acknowledgements](#references-license--acknowledgements)

---

## Quick Start

**1. Prerequisites** — install and build these first (details in [Dependencies & Environment Setup](#dependencies--environment-setup)):

- [Cloudy C25](https://gitlab.nublado.org/cloudy/cloudy) built as a static library (`libcloudy.a`)
- [HEASoft 6.33+](https://heasarc.gsfc.nasa.gov/docs/software/heasoft/) (provides Xspec model libraries)
- `g++` with C++17

**2. Edit the build paths in `Makefile`** — `CLOUDY_SRC` ships as a placeholder and `XSPEC_LIB` is hardcoded to the author's machine, so the build fails until you fix them (`CLOUDY_LIB` is auto-derived):

```make
CLOUDY_SRC = /path/to/cloudy/source
CLOUDY_LIB = $(CLOUDY_SRC)            # auto-derived; leave as is
XSPEC_LIB  = $(HEADAS)/lib            # set to your HEADAS (exported in step 4)
```

**3. Install the custom Cloudy continuum mesh** (required — see [resolution file](#resolution-file)):

```bash
cp /path/to/cloudy/data/continuum_mesh.ini /path/to/cloudy/data/continuum_mesh.ini.bak
cp resolution/smooth_hump.ini /path/to/cloudy/data/continuum_mesh.ini
```

**4. Set up the runtime environment and build:**

```bash
export HEADAS=/path/to/heasoft-6.33.2/<arch>
source $HEADAS/headas-init.sh        # must be sourced before every run
make
```

**5. Run a smoke test** (a `nthcomp`-illuminated slab):

```bash
./maindaocl -corona nthcomp -Gamma 2.0 -kT_e 60 -kT_bb 0.1 -nh 15 -zeta 3
```

Expected stdout (abridged) — when you see the parameter banner and a run hash, the run started correctly:

```
Corona:     nthcomp  Gamma=2.0000  kT_e=60.00 keV  kT_bb=0.1000 keV
Parameters: nh=15.0000  zeta=3.0000 (xi=1.0000e+03)  frac=-1.0000  ...
Run hash:   a1b2c3d4  →  results/a1b2c3d4/
```

Output lands in `results/<hash>/` (see [Outputs](#outputs)). **The first run is slow** because it builds and disk-caches the Compton kernel; later runs load it instantly.

---

## Dependencies & Environment Setup

| Software | Version | Purpose |
|----------|---------|---------|
| **Cloudy** | C25 | Atomic physics: emissivity, opacity, thermal & ionisation balance |
| **HEASoft / Xspec** | 6.33+ | Corona model libraries (`nthcomp`, `comptt`) |
| **g++** | C++17 | Compiler |
| **Python** | 3.8+ | Web UI (optional) |
| **Flask** | 2.0+ | Web UI server (optional: `pip install flask`) |

### Cloudy

Build Cloudy as a static library (`libcloudy.a`) and point `CLOUDY_SRC` / `CLOUDY_LIB` at its `source/` directory in the `Makefile`. DAO links against `gitlab.nublado.org/cloudy/cloudy` and calls it through its public API (`cdInit()` / `cdDrive()`).

<a name="resolution-file"></a>
**Continuum-mesh resolution file (required).** DAO ships a tailored Cloudy mesh at `resolution/smooth_hump.ini` that boosts the spectral resolving power around the iron-K region while staying modest elsewhere. Install it as Cloudy's `data/continuum_mesh.ini` (back up the original first, as in [Quick Start](#quick-start)). Without it, Cloudy's energy grid is too coarse to run the scattering kernel correctly. To choose a different resolution, edit the parameters at the top of `resolution/make_smooth_hump.py`, run it to regenerate `smooth_hump.ini`, then reinstall. See `resolution/README.md` for the parameters and file format.

### HEASoft / Xspec

Required for `nthcomp`, `comptt`, **and all production runs** (the executable links `-lXSFunctions`). Source it before every run, or dyld will fail to load the Xspec shared libraries:

```bash
export HEADAS=/path/to/heasoft-6.33.2/<arch>
source $HEADAS/headas-init.sh
```

Set `XSPEC_LIB = $(HEADAS)/lib` in the `Makefile` so the link-time and run-time libraries stay consistent.

### Compton kernel cache

The Compton kernel is precomputed once and cached on disk. The cache directory is set by `COMPTON_CACHE_DIR` (default: current working directory → `./kernel/`):

```bash
export COMPTON_CACHE_DIR=/path/to/shared/cache   # optional; kernels are large & parameter-independent
```

Files: `kernel/kernel_norm_NE{}_NI{}_NT{}.bin` (angle-dependent) and `kernel/avgkernel_norm_NE{}_NT{}.bin` (angle-mean, used by `-angsca false`). Grid dimensions are encoded in the filename so test and production grids coexist.

### Build

```bash
make           # build maindaocl
make clean     # remove object files and binaries
```

---

## Usage

The `-corona` flag is **required**, and each model requires specific parameters (see [required parameters](#required-parameters-per-corona-model)).

```bash
# Cutoff power law
./maindaocl -corona cutoffpl -Gamma 2.0 -Ecut 300 -nh 16 -zeta 4 -frac 0.5

# Thermal Comptonisation (nthcomp)
./maindaocl -corona nthcomp -Gamma 2.0 -kT_e 60 -kT_bb 0.1 -nh 15 -zeta 3

# CompTT, with enhanced iron
./maindaocl -corona comptt -kT_e 50 -kT_bb 0.05 -taup 1.0 -Afe 3.0

# Test mode — compPS benchmark slab (see below)
./maindaocl -test_rt compps -corona blackbody -kT_e 60 -kT_bb 0.1 -tau 0.5
```

### Parameter Reference

All `kT_*` temperatures are in **keV**. Defaults are from `source/params.cpp`.

| Flag | Default | Description |
|------|---------|-------------|
| `-corona` | **required** | Corona model: `powerlaw`, `cutoffpl`, `nthcomp`, `comptt`, `blackbody` |
| `-Gamma` | unset | Photon index Γ (powerlaw, cutoffpl, nthcomp) |
| `-Ecut` | unset | High-energy cutoff, keV (cutoffpl) |
| `-E_low_cut` | 0.1 | Low-energy exponential cutoff, keV (powerlaw, cutoffpl) |
| `-kT_e` | unset | Electron temperature, keV (nthcomp, comptt; also slab Tₑ in compps test mode) |
| `-kT_bb` | unset | Seed-photon temperature, keV (nthcomp, comptt, blackbody) |
| `-taup` | unset | Plasma optical depth (comptt) |
| `-nh` | 15 | log₁₀ hydrogen density [cm⁻³] |
| `-zeta` | 3.0 | Ionisation parameter: ξ = 10^zeta (see convention below) |
| `-frac` | **−1** | Flux ratio F_corona / F_disk; **≤0 (default) → corona only, no disk component** |
| `-incidence` | 0.7071 | cos θ of corona incidence (snapped to nearest Gauss–Legendre node) |
| `-kT_disk` | 0.35 | Disk blackbody temperature, keV (thermal component illuminating from below) |
| `-Afe` | 1.0 | Iron abundance [solar] |
| `-angsca` | true | `true`/`1`/`yes` → angle-dependent kernel; `false`/`0`/`no` → angle-mean kernel |
| `-test_rt <mode>` | off | Test mode (skips Cloudy). Mode token **required**; only public mode is `compps` |
| `-tau` | 0.5 | Slab vertical Thomson optical depth (compps test mode) |

**Ionisation parameter convention.** DAO uses **ξ = (4π)² J / n_H**, i.e. the mean intensity is normalised as **J = ξ·n_H / (4π)²**, with `zeta = log₁₀ ξ`. (This is the convention implemented in `source/radiation.cpp`.)

<a name="required-parameters-per-corona-model"></a>
**Required parameters per corona model:**

| Model | Requires |
|-------|----------|
| `powerlaw` | `-Gamma` |
| `cutoffpl` | `-Gamma -Ecut` |
| `nthcomp` | `-Gamma -kT_e -kT_bb` |
| `comptt` | `-kT_e -kT_bb -taup` |
| `blackbody` | `-kT_bb` |

### Test mode: compPS benchmark

`-test_rt compps` (the command in the [Usage examples](#usage) above) runs an **isothermal, pure-scattering** slab seeded by a **bottom blackbody**, benchmarked against Xspec's compPS (Poutanen & Svensson 1996). The slab temperature is `-kT_e` [keV] and its vertical Thomson depth is `-tau`. This mode skips Cloudy entirely and uses a double-Gauss angle grid; it is a validation benchmark, not a production reflection run.

---

## Outputs

Results are written to `results/<hash>/`, where `<hash>` is an 8-character FNV-1a hash of the physics parameters. **The hash is deterministic** — identical parameters reuse (and overwrite) the same directory. Every run directory contains a `params.json` with the full parameter set.

| File | When | Contents |
|------|------|----------|
| `params.json` | always | Full parameter set (JSON) |
| `emergent_iter{NNN}.dat` | production | Emergent surface intensity: E, I_corona, I_disk, I(μ) per angle |
| `moments_iter{NNN}.dat` | production | Angular moments J0, J2, J3 at all depths and energies |
| `profile_iter{NNN}.dat` | production | Depth profile: τ, T, n_e, heating, cooling, log ξ |
| `emergent_compps.dat` | compps test | Emergent spectrum (no `_iter` files) |

### Data units (per-eV throughout)

| Quantity | Units |
|----------|-------|
| Energy E | eV |
| j_ν (emissivity) | erg cm⁻³ s⁻¹ eV⁻¹ |
| κ_abs, κ_sct (opacity) | cm⁻¹ |
| I_ν (specific intensity) | erg cm⁻² s⁻¹ eV⁻¹ sr⁻¹ |
| I_corona, I_disk | erg cm⁻² s⁻¹ eV⁻¹ |

---

## Web UI

A browser-based configurator and results viewer (requires `pip install flask`):

```bash
python3 ui.py        # serves http://127.0.0.1:5200
```

| Page | URL | Description |
|------|-----|-------------|
| Configurator | `/` | Parameter sliders, corona-model selector, command generator |
| Results viewer | `/plots` | Emergent spectra, mean intensity, temperature profiles |
| Convergence | `/convergence` | Temperature-convergence history across outer iterations |
| Parameter reference | `/docs` | Full parameter documentation, references & license |

The configurator generates commands to copy into your terminal; the results viewer reads `results/*/params.json` and the `.dat` files.

---

## How It Works

DAO solves the angle- and energy-dependent transfer equation in a plane-parallel disc atmosphere, illuminated by an X-ray corona from above and a thermal disc from below. At each depth point, **Cloudy** returns the emissivity, absorption opacity, and thermal/ionisation structure; the **RT solver** then propagates the radiation field through the slab using a **cubic Bézier short-characteristics formal solution with Lambda iteration**, including **exact, angle/energy-dependent Compton redistribution**. The two are iterated to convergence.

### Execution flow

```
CLI arguments
   │
   ▼
read_params()  →  ModelParams (validate corona model + required params)
   │
   ▼
init_rt_grids()  →  angle (Gauss–Legendre), depth (Thomson τ), energy grids
   │
   ▼
rad.illum.compute()  →  corona spectrum (dispatch) + disk blackbody, normalised by ξ
   │
   ▼
Precompute Compton kernel  (load from disk cache or compute + save)
   │
   ▼
Dispatch:
   ├── test mode  →  run_test_rt()  (compps benchmark slab; synthetic
   │                 1000-bin log grid 0.01–1000 keV, no Cloudy)
   │
   └── production →  Bootstrap Cloudy → adopt its energy mesh (~3300 RT bins)
                     Outer loop:
                       ├── Cloudy depth sweep → j_ν, κ_abs, κ_sct
                       ├── RT solve (Bézier formal solution + Lambda iteration)
                       ├── Compute moments J and ionisation parameter ξ
                       ├── Check convergence (max |ΔT/T|, max |Δ log ξ|)
                       └── Save results
```

### Corona-model dispatch

| Model | Spectral shape | Interface |
|-------|----------------|-----------|
| `powerlaw` | E^(1−Γ)·exp(−E_lo/E) | analytic |
| `cutoffpl` | E^(1−Γ)·exp(−E/E_cut)·exp(−E_lo/E) | analytic |
| `nthcomp` | thermal Comptonisation | Xspec `donthcomp_()` |
| `comptt` | Comptonisation | Xspec `C_compTT()` |
| `blackbody` | Planck function | analytic |

### Code structure

```
DAO/
├── maindaocl.cpp          Main entry point
├── Makefile               Build configuration
├── ui.py                  Web UI (Flask)
├── source/                C++ modules:
│   ├── compton_kernel.*           Exact azimuth-integrated kernel (Madej+ 2017)
│   ├── avg_compton_kernel.*       Angle-mean kernel (-angsca false)
│   ├── compton_cross_section.*    Relativistic σ(E,T) (Poutanen & Svensson 1996)
│   ├── compton_rt.*               RT solver (Bézier short characteristics + Λ iteration)
│   ├── source.*                   Source function: thermal + Compton scattering
│   ├── radiation.*                Radiation arrays, moments, convergence
│   ├── corona_models.*            Corona spectrum dispatch
│   ├── cloudy_interface_v2.cpp    Cloudy emissivity/opacity extraction
│   ├── rt_grids.*                 Angle / depth / energy grids
│   ├── params.*                   CLI parsing, validation, run hash
│   ├── save_results.*             Output to results/<hash>/
│   ├── production.*               Cloudy ↔ RT outer iteration
│   └── test_rt.*                  compps benchmark mode
├── resolution/            Custom Cloudy continuum-mesh files + generator
├── paper_data/            Data for the DAO paper
├── image/                 UI assets
├── kernel/                Cached Compton kernels (binary, generated at runtime)
└── results/              Run outputs: results/<hash>/
```

---

## References, License & Acknowledgements

### References

If you use DAO in published work, please cite the **DAO paper**, **Madej et al. (2017)** for the Compton kernel, and **Gunasekera et al. (2025)** for Cloudy. A [`CITATION.cff`](CITATION.cff) is provided.

- **Atomic physics (Cloudy):** Gunasekera et al. 2025, arXiv:2508.01102 (Cloudy C25).
- **Compton kernel:** Madej, Różańska, Majczyna & Należyta 2017, MNRAS 469, 2032; Nagirner & Poutanen 1993, A&A 275, 325.
- **Cross section & A23 normalisation:** Poutanen & Svensson 1996, ApJ 470, 249.
- **RT solver:** Auer 2003, ASP Conf. Ser. 288, 3; de la Cruz Rodríguez & Piskunov 2013, ApJ 764, 33; Hubeny & Mihalas 2015, *Theory of Stellar Atmospheres* §12.4; Suleimanov, Poutanen & Werner 2012, A&A 545, A120.
- **Corona models:** Zdziarski, Johnson & Magdziarz 1996, MNRAS 283, 193 (nthcomp); Titarchuk 1994, ApJ 434, 570 (comptt); Arnaud 1996, ASP Conf. Ser. 101, 17 + HEASoft (XSPEC).
- **Ionisation parameter:** Tarter, Tucker & Salpeter 1969, ApJ 156, 943.

### License

DAO's source code is released under the **MIT License** (see [`LICENSE`](LICENSE)). DAO links an **unmodified**, separately installed Cloudy through its public API and ships no Cloudy source. A few physics steps are reimplemented on DAO's own grid: per-line escape probabilities (calling Cloudy's `rt_escprob`), subtraction of Cloudy's bound-electron Compton-recoil opacity (reproducing `opacity_addtotal.cpp`), and inner-shell fluorescence (reproducing `prt_lines.cpp`).

Third-party dependencies, obtained independently and each under its own license:

| Component | License | Source |
|-----------|---------|--------|
| **Cloudy** (G. J. Ferland et al.) | [zlib](https://opensource.org/licenses/Zlib) | <https://gitlab.nublado.org/cloudy/cloudy> |
| **HEASoft / XSPEC** (NASA HEASARC) | NASA open-source (HEASARC) | <https://heasarc.gsfc.nasa.gov/docs/software/heasoft/> |
| **Compton redistribution kernel** | C++ port of J. Madej's public Fortran | <https://www.astrouw.edu.pl/~jm/software.html> |

### Acknowledgements

We gratefully acknowledge **Prof. Jerzy Madej** (University of Warsaw) for his publicly available Fortran implementation of the exact Compton redistribution function, on which our C++ port is based.

For questions, suggestions, or collaboration, open an issue on GitHub or contact the author at huangym23@m.fudan.edu.cn.
