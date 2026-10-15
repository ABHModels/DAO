# DAO — X-ray Reflection Spectroscopy Model

**Version 1.0.0**

<img align="right" src="image/bg_blackhole.png" alt="Black hole and accretion disk" width="300">

DAO computes rest-frame X-ray reflection spectra using Cloudy for atomic
physics, the XSPEC model library for incident continua, and angle-dependent
Compton radiative transfer. Gas temperatures are determined from radiative
energy balance.

**Author:** Yimin Huang · Fudan University; University of Bristol · yi-min.huang@bristol.ac.uk

<br clear="all">

## 1. Install the dependencies

Use macOS or Linux with a C++17 compiler (`g++` or `clang++`), `make`, and Git.
Install these two packages first:

- **[Cloudy C25](https://gitlab.nublado.org/cloudy/cloudy):** download the source
  and its required atomic data following the upstream instructions. Build it
  in `source/`; DAO needs `libcloudy.a`, the headers, and `cloudyconfig.h`.
- **[HEASoft / XSPEC](https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/install.html):**
  include XSPEC when installing. DAO links to its model library, `libXSFunctions`,
  for `nthcomp` and `compTT`. The library is required to link the DAO executable.

The commands below use bash/zsh. Replace the example paths with your own:

```bash
export DAO_CLOUDY_ROOT=/absolute/path/to/cloudy
export HEADAS=/absolute/path/to/heasoft/platform-directory
source "$HEADAS/headas-init.sh"

make -C "$DAO_CLOUDY_ROOT/source" -j4
```

Use compatible C++ compilers for Cloudy and DAO. `HEADAS` is the installed
platform directory containing `headas-init.sh` and `lib/`, not the top-level
HEASoft source directory. Initialize HEASoft in every new terminal used for DAO.

## 2. Download DAO

```bash
git clone https://github.com/ABHModels/DAO.git
cd DAO
```

Run the remaining commands from this directory.

## 3. Configure Cloudy's runtime files — required

DAO requires the supplied [`config/xray12.ini`](config/xray12.ini). It selects
the 12 elements used by DAO: H, He, C, N, O, Ne, Mg, Si, S, Ar, Ca, and Fe.
This file and the continuum mesh are supplied in `config/`.

**DAO takes its production energy grid from Cloudy.** The supplied default,
[`config/continuum_mesh.ini`](config/continuum_mesh.ini), uses logarithmic
spacing over 1 eV–1000 keV. You can replace it with your own Cloudy-compatible
mesh, keeping the filename **`continuum_mesh.ini`**. From the DAO directory, run:

```bash
export CLOUDY_DATA_PATH="$PWD/config:.:$DAO_CLOUDY_ROOT/data"
```

**The mesh must be named `continuum_mesh.ini`.**

## 4. Compile DAO

Override the paths supplied in the Makefile with your installation paths:

```bash
make -j4 CLOUDY_SRC="$DAO_CLOUDY_ROOT/source" \
         CLOUDY_LIB="$DAO_CLOUDY_ROOT/source" \
         XSPEC_LIB="$HEADAS/lib"
```

This creates `./maindaocl`. If Cloudy was built with `clang++`, add
`CXX=clang++` to this command. Alternatively, save these paths in the Makefile.
After changing compilers or dependency installations, run `make clean` before
rebuilding. A missing `libcloudy.a` means Cloudy has not been built at the given
path; a missing XSPEC library usually means `HEADAS` or `XSPEC_LIB` is incorrect.

> **Compiler compatibility:** Build Cloudy and DAO with the same C++ compiler
> and version. The author's tested environment is Apple Clang 21.0.0 on macOS,
> using C++17; other compiler environments have not yet been validated.
> If you encounter compilation problems, please contact the author with your
> compiler version, build command, and error output.

## 5. Run a reflection model

Set a dedicated kernel-cache directory (recommended: `DAO_KERNEL_TABLE`):

```bash
mkdir -p "$PWD/DAO_KERNEL_TABLE"
export COMPTON_CACHE_DIR="$PWD/DAO_KERNEL_TABLE"
export DAO_KERNEL_THREADS=4
```

Environment variables used for setup and runtime:

| Variable | Default when unset | Meaning |
|----------|--------------------|---------|
| `DAO_CLOUDY_ROOT` | None; set for the setup commands | Cloudy installation root, containing `source/` and `data/`. |
| `HEADAS` | None; set when initializing HEASoft | HEASoft platform directory, containing `headas-init.sh` and `lib/`. |
| `CLOUDY_DATA_PATH` | Cloudy's build-configured search path | Set to `$PWD/config:.:$DAO_CLOUDY_ROOT/data`; keep `.` so Cloudy can read DAO's temporary incident spectra. |
| `COMPTON_CACHE_DIR` | Current working directory (`.`) | Directory for loading and saving Compton kernels. Create it first; `DAO_KERNEL_TABLE` is a recommended name. |
| `DAO_CLOUDY_WORKERS` | Up to 4, depending on available CPU threads | Parallel Cloudy processes, 1–64. Lower this if memory is limited. |
| `DAO_RT_THREADS` | Up to 16, depending on available CPU threads | Threads for parallel radiation calculations, 1–256. |
| `DAO_KERNEL_THREADS` | Up to 16 for kernel construction; up to 8 for normalization, depending on available CPU threads | Threads for kernel generation, 1–256. Setting this does not speed up loading an existing cache. |

CPU defaults use the logical thread count reported by the system; set worker
counts explicitly to match your cluster job's CPU allocation. Use an absolute
`COMPTON_CACHE_DIR` to share a cache between runs. Without this variable, cache
files are saved directly in `.`; DAO does not automatically create `./kernel/`.
`resume_model.py` also checks for `DAO_KERNEL_TABLE` beside the checkout and uses
it if present, without overriding an explicitly set `COMPTON_CACHE_DIR`.

```bash
DAO_CLOUDY_WORKERS=4 DAO_KERNEL_THREADS=4 DAO_RT_THREADS=4 \
./maindaocl -corona nthcomp -Gamma 2.0 -kT_e 60 -kT_bb 0.01 \
  -nh 15 -zeta 3 -frac -1 -incidence 0.7071 -angsca 1
```

This illuminates a constant-density slab with `nthcomp`, using
nH = 10¹⁵ cm⁻³, log ξ = 3, and angle-dependent Compton scattering.
The default slab has 48 logarithmic depth cells and reference Thomson depth 5.
The first run builds a Compton-kernel cache and can take substantially longer
than subsequent runs. These are full atmosphere calculations, not quick tests.
Reduce the worker counts to `1` if memory is limited.

All command-line options accepted by `./maindaocl`:

| Option | Default | Meaning |
|--------|---------|---------|
| `-corona MODEL` | Required | Incident spectrum: powerlaw, cutoffpl, nthcomp, comptt, or blackbody. |
| `-Gamma G` | Required when used | Power-law slope for powerlaw, cutoffpl, or nthcomp; larger values give a softer spectrum. |
| `-Ecut E` | Required for cutoffpl | High-energy cutoff in keV. |
| `-E_low_cut E` | 0.1 | Low-energy cutoff in keV, for powerlaw and cutoffpl. |
| `-kT_e T` | Required when used | Electron temperature in keV for nthcomp/comptt; slab temperature in RT tests. |
| `-kT_bb T` | Required when used | Seed-photon temperature in keV for nthcomp/comptt, or temperature of the incident blackbody. |
| `-taup T` | Required for comptt | Optical depth of the Comptonizing corona. |
| `-nh N` | 15 | Log₁₀ hydrogen density in cm⁻³: 15 means 10¹⁵ cm⁻³. |
| `-zeta Z` | 3 | Log₁₀ ionization parameter ξ: 3 means ξ = 1000 erg cm s⁻¹. |
| `-frac R` | -1 | Disk/corona flux ratio in the normalization band; R ≤ 0 gives corona-only illumination. |
| `-kT_disk T` | 0.35 | Bottom blackbody temperature in keV; used when -frac is positive. |
| `-incidence MU` | 0.7071 | Cosine of the angle from the slab normal; snapped to the nearest angle node. Use -2 for all downward directions. |
| `-angsca 0/1` | 1 | 1 keeps the scattering-angle dependence; 0 uses the angle-averaged kernel. |
| `-Afe A` | 1 | Iron abundance relative to solar: 1 is solar, 2 is twice solar. |
| `-O [0/1]` | Off | Save final oxygen ion fractions; -O alone enables this output. |
| `-Fe [0/1]` | Off | Save final iron ion fractions; -Fe alone enables this output. |
| `-max_outer_iter N` | 100 | Total iteration limit, including any completed iterations loaded for continuation. |
| `-resume_iter N` | Off | Continue from the last completed iteration N with the same model settings; resume_model.py supplies these settings for you. |
| `-resume_check` | Off | Check a saved state and exit without calculating; requires -resume_iter. |
| `-test_rt MODE` | Off | Run a scattering-only test: compps uses bottom blackbody illumination; test_avg uses top coronal illumination. |
| `-tau T` | 0.5 | Vertical Thomson depth for the compps test; does not change the production slab depth. |
| `-label "TEXT"` | Empty | A name shown with the run; does not change the physics. |
| `-verbose [0/1]` | Off | Include detailed diagnostics in the terminal and run.log; -verbose alone enables them. |

Options that take 0/1 also accept false/true or no/yes.
Required inputs depend on the selected spectrum; unused spectrum parameters can be omitted.

## 6. Find the results

The startup message prints `results/<hash>/`. The terminal then shows one
summary per completed outer iteration. The same summaries go to `run.log`;
detailed diagnostics are included only with `-verbose`.

| File | Contents |
|------|----------|
| `params.json`, `RUN.txt` | Input parameters and run information |
| `thermal_status.json` | Run state; **`converged`** identifies a successful solution |
| `emergent_iterNNN.dat` | Upper-surface spectra at the outgoing angular-grid nodes |
| `profile_iterNNN.dat` | Depth profiles, including temperature and electron density |
| `thermal_history.dat` | Convergence history and outgoing/incoming flux ratio |
| `thermal_budget_iterN.dat` | Local and whole-slab energy budgets |
| `O/ion_fractions.dat`, `Iron/ion_fractions.dat` | Optional ion fractions, saved after convergence |

To save oxygen and iron ion fractions versus depth, append `-O 1 -Fe 1` to the
model command. Either flag can be enabled independently; both default to off.
DAO reads the ion populations directly from Cloudy and creates
`results/<hash>/O/ion_fractions.dat` and
`results/<hash>/Iron/ion_fractions.dat` only after convergence. The tables include
all ion stages, from neutral to fully stripped, together with depth, temperature,
and electron density.

Intermediate spectra can exist even if a run fails to converge. Use the final
iteration identified by `thermal_status.json`. Both slab faces are included in
the reported energy balance. Repeating identical physical parameters reuses
the same result directory and can overwrite its outputs.

Some models, especially at very low ionization, need more than 100 iterations.
Although `-max_outer_iter` sets the limit, the iterations needed for convergence
are not known in advance. `resume_model.py` continues from saved results without
starting over. Replace `YOUR_HASH` with the model's run hash:

```bash
python3 resume_model.py results/YOUR_HASH --max-iterations 300
```

Here 300 is the total iteration limit, including completed iterations;
the run stops earlier if it converges. Add `--check` to validate the checkpoint
without running the model.

Options for `resume_model.py`:

| Option | Default | Meaning |
|--------|---------|---------|
| `results/YOUR_HASH` | Required | Existing run directory to continue. |
| `--max-iterations N` | 300 | Total iteration limit, including completed iterations. |
| `--iteration N` | Last completed step | Choose the saved step; it must be the last completed iteration. |
| `--check` | Off | Validate the saved state without running the model. |
| `--dry-run` | Off | Print the continuation command without running it. |
| `-h / --help` | — | Show usage and exit. |

For a browser-based launcher and spectrum viewer, start this from the same
initialized terminal:

```bash
python3 -m pip install flask numpy matplotlib
python3 ui.py
```

Open **http://127.0.0.1:5200**. Python is optional for command-line model runs.

## Known Issues
When ionization is very low, the thermal equilibirum is unstable at lowest depth. A fixed version is testing, we will publish when it is stable.

## Citation and license

Please cite the DAO paper and the underlying methods, particularly
[Madej et al. (2017)](https://ui.adsabs.harvard.edu/abs/2017MNRAS.469.2032M),
[Cloudy C25](https://doi.org/10.48550/arXiv.2508.01102), and
[XSPEC](https://ui.adsabs.harvard.edu/abs/1996ASPC..101...17A).
See [CITATION.cff](CITATION.cff) for citation metadata.

DAO's original code uses the **MIT License**; see [LICENSE](LICENSE).
Cloudy is used unmodified under its zlib license. The separately installed
XSPEC model library retains its upstream terms. The C++ adaptations of
Jerzy Madej's Compton routines are included with his permission as a co-author
of the DAO paper. Neither the Cloudy nor XSPEC library is bundled with DAO.
