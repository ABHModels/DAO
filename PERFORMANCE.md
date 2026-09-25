# Performance validation

This change preserves the kernel format, double precision, quadrature, band
thresholds, per-cell summation order, grids, and RT convergence criteria.
Cloudy sources and the linked Cloudy library were not modified.

## Implementation

- Compute independent kernel rows concurrently; retain exact row values instead
  of evaluating the retained band twice.
- Normalize one temperature at a time through a private mapping of a temporary
  row spool. Write completed temperatures sequentially and release private pages.
- Publish the complete cache by atomic rename and retain a demand-paged mapping.
- Parallelize independent normalization integrals and RT source-function depths.
- Reuse identical detailed-balance exponentials across angles without changing
  subsequent arithmetic or reduction order.

`DAO_KERNEL_THREADS` controls construction and normalization (default up to 8).
`DAO_RT_THREADS` controls source-function evaluation (default up to 16).
Both accept 1–256. More threads are not always faster.

## Recorded local evidence

Measured on macOS, Apple M5 Pro, 24 GB RAM, C++17 `-O3`, without fast-math.
Results are workload-specific, not universal speed or memory guarantees.

The final memory implementation was compared with the already accelerated
implementation before temperature-streamed normalization:

| Measurement | Before memory optimization | Final |
| --- | --- | --- |
| 2048-energy, 8-angle, 3-temperature cold construction, median of 3 | 8.426 s | 6.444 s |
| Observed peak process RSS in those runs | 354.84–354.86 MiB | 181.81–186.69 MiB |
| RSS immediately after initialization | about 354.85 MiB | 3.56–3.58 MiB |
| Full cold compPS run, 4 alternating pairs, median | 6.890 s | 6.615 s |

Initialization RSS is not full RT memory: mapped pages become resident when
accessed. OS file-cache memory and temporary disk space must also be considered.

Binary caches and emergent spectra compared byte-for-byte (zero tolerance).
Directional and angle-mean compPS runs also had identical convergence sequences.
Earlier RT validation compared complete in-memory intensity/mean-intensity
arrays for mixed temperatures, absorption/emission, and nonuniform depths at
1/3/15 threads; source arrays were compared at 1/2/8/15 threads.
Storage/worker assertions passed AddressSanitizer and UndefinedBehaviorSanitizer.

## Reproduction checks

```sh
make test_kernel_storage
./test_kernel_storage
make -j4 CLOUDY_SRC=/path/to/cloudy/source XSPEC_LIB=/path/to/heasoft/lib
```

Run each version in a separate empty run/cache directory; do not overwrite
existing spectra or reuse caches when measuring cold construction:

```sh
COMPTON_CACHE_DIR="$PWD" /path/to/maindaocl \
  -test_rt compps -corona blackbody -kT_e 60 -kT_bb 0.1 \
  -tau 0.5 -nh 12 -zeta 3 -frac 100 -incidence 0.5 -angsca 1
```

Repeat with `-angsca 0`. Configure runtime library search paths for the local
HEASoft installation. Compare corresponding caches and spectra with `cmp`.
The storage test alone is not a physics-accuracy certification.

The full 50-temperature production kernel, full coupled Cloudy calculation,
and cross-platform regression were not rerun. Byte identity of these tested
cases does not establish the physical accuracy of the original approximations.
