# DAO — X-ray Reflection Spectroscopy Model

**Version 1.0.0**

DAO computes rest-frame X-ray reflection spectra using Cloudy for atomic
physics, the XSPEC model library for incident continua, and angle-dependent
Compton radiative transfer. Gas temperatures are determined from radiative
energy balance.

**Author:** Yimin Huang · Fudan University; University of Bristol · yi-min.huang@bristol.ac.uk

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
This file and the continuum mesh are supplied in `config/`; there is no need
to copy them into Cloudy's installation.

**DAO takes its production energy grid from Cloudy.** The supplied default,
[`config/continuum_mesh.ini`](config/continuum_mesh.ini), uses logarithmic
spacing over 1 eV–1000 keV. You can replace it with your own Cloudy-compatible
mesh, keeping the filename **`continuum_mesh.ini`**. From the DAO directory, run:

```bash
export CLOUDY_DATA_PATH="$PWD/config:.:$DAO_CLOUDY_ROOT/data"
```

Cloudy searches these directories in order: `config/` supplies DAO's mesh and
initialization file, `.` lets it read the temporary incident spectra generated
in the current directory, and its installed `data` directory supplies the
remaining atomic data. The same configuration is used by DAO's parallel
Cloudy workers. Keep `.` in the search path so those spectra can be read.

The mesh must be named **`continuum_mesh.ini`** and its directory must come first
in `CLOUDY_DATA_PATH`. To use a custom mesh, replace `config/continuum_mesh.ini`
with your own before starting the model. The mesh is
read at runtime; changing it does not require recompiling Cloudy.

Keep the initialization filename **`xray12.ini`**: DAO loads it with Cloudy's
`init "xray12.ini"` command. This file is required even when using a custom mesh.

Keep the selected configuration files unchanged during a run. The environment
variable applies to programs started from this terminal; existing processes
retain their own environment. In each new terminal, set `DAO_CLOUDY_ROOT`,
initialize HEASoft, and export `CLOUDY_DATA_PATH` again from the DAO directory.
See [resolution/README.md](resolution/README.md) for custom meshes.

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

## 5. Run a reflection model

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

| Option | Meaning |
|--------|---------|
| `-nh 15` | log₁₀ hydrogen number density in cm⁻³, not column density |
| `-zeta 3` | log₁₀ ionization parameter ξ |
| `-Gamma 2` | Incident photon index |
| `-kT_e 60`, `-kT_bb 0.01` | Coronal electron and seed-photon temperatures in keV |
| `-frac -1` | Corona-only illumination; no lower-boundary disk component |
| `-incidence 0.7071` | Incident cosine, mapped to the nearest angular-grid node |
| `-incidence -2` | Isotropic illumination over the downward hemisphere |
| `-angsca 1` / `0` | Angle-dependent / angle-averaged Compton scattering |
| `-Afe 1` | Iron abundance relative to solar |
| `-O 1 -Fe 1` | Save final oxygen and iron ion fractions |
| `-verbose 1` | Also show detailed diagnostics in the terminal |

For corona-only illumination, the flux integrated over the DAO energy grid is ξ nH / (4π).

## 6. Find the results

The startup message prints `results/<hash>/`. The terminal then shows one
summary per completed outer iteration; detailed diagnostics go to `run.log`.

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
**No additional environment variables or Cloudy output-path settings are needed.**
DAO reads the ion populations directly from Cloudy and creates
`results/<hash>/O/ion_fractions.dat` and
`results/<hash>/Iron/ion_fractions.dat` only after convergence. The tables include
all ion stages, from neutral to fully stripped, together with depth, temperature,
and electron density.

Intermediate spectra can exist even if a run fails to converge. Use the final
iteration identified by `thermal_status.json`. Both slab faces are included in
the reported energy balance. Repeating identical physical parameters reuses
the same result directory and can overwrite its outputs.

For a browser-based launcher and spectrum viewer, start this from the same
initialized terminal:

```bash
python3 -m pip install flask numpy matplotlib
python3 ui.py
```

Open **http://127.0.0.1:5200**. Python is optional for command-line model runs.

## Further reading

- [Paper figures, benchmark data, and plotting instructions](paper_data/README.md)
- [Release notes](CHANGELOG.md)

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
