# Low-temperature Compton-kernel validation

Validated on 2026-09-26 against public `main` at `a05d541`. The change keeps
the exact electron redistribution profile, Klein-Nishina cross section, 32-point
Gauss-Laguerre electron integral, and `1.21*nH` scattering density. The
approximate directional mode now uses its selected approximate profile, which
the previous angular integrator ignored. It changes
the angular and energy-cell quadrature and the discrete normalization. Its
numerical outputs therefore are **not bitwise identical** to the previous
kernel; the old values can be inaccurate where a narrow peak falls between
quadrature nodes. No claim of an identical emergent spectrum is made here.

## Independent angular reference

`source/kernel_reference_probe.cpp` evaluates the same physical profile with
256- and 512-point Gauss-Legendre angle rules. It requires their relative
difference to be below `1e-4` and the production kernel to agree with the
512-point result within `2e-3`, for both directional and angle-mean kernels.
It covers exact (`ktype=1`) and approximate (`ktype=0`) profiles. The largest
measured production error across seven cases was `1.452e-5`; the largest
256-versus-512 reference difference was `7.289e-8`.

| Case | Temperature (K) | Production directional | Dense directional | Production angle mean | Dense angle mean | Largest relative error |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| hot low energy | 1e9 | 59.0442974659 | 59.0442978155 | 74.9786264718 | 74.9775375309 | 1.452e-5 |
| hot high energy | 1e9 | 1.37609385748 | 1.37609385748 | 1.79257971451 | 1.79257971424 | 1.535e-10 |
| hot cross angle | 1e9 | 0.328125269657 | 0.328125269657 | 1.09110151434 | 1.09110151435 | 1.004e-11 |
| cold recoil | 1e4 | 4.26357462355 | 4.26357462356 | 2.11851689112 | 2.11851689112 | 2.913e-12 |
| warm upscatter | 1e6 | 0.397280972092 | 0.397280972092 | 0.0668170905969 | 0.0668170905969 | 1.662e-13 |
| cold approximate | 1e4 | 0.00149034845124 | 0.00149034845121 | 0.0194536055550 | 0.0194536055550 | 2.149e-11 |
| warm approximate | 1e6 | 37.1447417573 | 37.1447417573 | 49.8493770325 | 49.8493770325 | 1.136e-13 |

The earlier public kernel returned `68.3654050773` instead of `59.0442978155`
for the hot low-energy directional case (15.8% high), `1.64288436313`
instead of `1.37609385748` for the hot high-energy directional case (19.4%
high), and `355.556718370` instead of `4.26357462356` for the cold recoil
directional case. These differences are why equality to the old discrete
values is not an appropriate physics criterion.

For a direct baseline report, build this probe against the source files at
`a05d541` and run it with `--report-only`; that flag prints discrepancies
without applying the new-kernel acceptance threshold. The dense reference
uses the probe and quadrature header from this change in both builds.

## Discrete redistribution and source function

`source/test_kernel_low_temp.cpp` passed with the 32-point electron integral.
It checks finite and positive kernels, directional and angle-mean cross-section
sum rules to `2e-3`, directional detailed balance to `1e-10`, finite source
functions in both modes, a cold single-bin row, smooth high-energy source
values on a 24-bin mesh, and cache reload/grid-mismatch handling. The smoothness
criterion is that adjacent tested source values do not rise by more than 3%.
This is a kernel-level and one-source-evaluation regression, not a full
radiative-transfer spectrum comparison.

`source/test_kernel_storage.cpp` also passed: ordered row spooling, private
mapping, atomic replacement, truncated payload rejection, and worker cleanup.
`maindaocl` compiled with the local Cloudy and HEASoft libraries. The compiler
reported one pre-existing unused-function warning in
`source/cloudy_interface_v2.cpp`.

Reproduce the standalone checks with:

```bash
make test_kernel_reference test_kernel_low_temp test_kernel_storage
./test_kernel_reference
./test_kernel_low_temp
./test_kernel_storage
```

The full Cloudy iteration and emergent spectrum were not run for this
validation. They remain the final end-to-end comparison before a release.
