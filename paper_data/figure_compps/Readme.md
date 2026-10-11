# Comptonised slab benchmark against compPS

The current reproduction is `compps_compare_a67d34d8.pdf`
(and the matching PNG). Run `a67d34d8` uses the same
`compton_rt_solve()` as production and `test_avg`: constant-cell analytic
transport, cell-volume-averaged intensities, photon-conserving scattering
operators, and a source-iteration tolerance of 1e-7 for three consecutive steps.

## Parameters and angular grid

- kTe = 60 keV, kTbb = 0.1 keV, vertical Thomson depth = 0.5.
- Isothermal pure scattering, blackbody injection from the bottom.
- 48 logarithmic depth cells, 1000 logarithmic energies over 0.01--1000 keV.
- Ten double-Gauss angular points, matching the original paper figure.
  This benchmark build changes only angular resolution; the default production
  executable still has eight angular points and uses exactly the same solver.
- compPS spectra are the original curated reference files, copied byte-for-byte
  after checking kTe, kTbb, tau and all angular nodes. SHA256 checksums and original
  paths are recorded in `reference_provenance.json`.
- The extreme cosines 0.046910 and 0.953090 are clamped to 0.05 and 0.95 by
  compPS, as in the original figure. The three interior nodes match exactly.

This directory contains only the current benchmark dataset. The previous
DAO spectra and figure have been moved out of the active dataset and preserved
in a verified local backup. The script reads the local `params.json` and
`emergent_compps.dat` by default.

## Reproduce the DAO calculation

From the repository root, with HEASoft libraries on the runtime library path:

```bash
make -j4 maindaocl compps_native_benchmark
DAO_KERNEL_THREADS=4 DAO_RT_THREADS=4 ./results/compps_native/maindaocl \
  -test_rt compps -corona blackbody -kT_e 60 -kT_bb 0.1 -tau 0.5 \
  -nh 12 -zeta 3 -frac 100 -incidence 0.5 -angsca 1
```

The alternate executable is compiled from the same sources using
`-DDAO_RT_ANGLES=10`. No alternate transfer algorithm is compiled or selectable.
New test hashes identify both the shared solver and angular resolution, so they
cannot overwrite the former Bezier results. Test modes prescribe temperature;
they do not run the production temperature search.

## Reproduce the figure

```bash
MPLBACKEND=Agg python paper_data/figure_compps/comppsGenerator.py --reuse-reference
```

This uses the saved DAO spectrum and original compPS reference files without
requiring PyXspec. To regenerate compPS instead, initialize HEASoft and omit
`--reuse-reference`. Add `--show` for an interactive window. `--data-dir` and
`--output-dir` select other datasets/destinations. The original layout, colours,
line styles and normalization are retained: each left-panel spectrum is divided
by its own integrated flux; the right panel plots the limb law using mu*I.

## Validation

The shared solver passes analytic absorption/emission, photon conservation,
Wien equilibrium, failure-on-nonconvergence, both angular modes, and a four-cell
Cloudy integration test. For the full ten-angle benchmark, the relative photon
number error is -5.9841e-14 and the signed boundary/volume energy identity error
is -2.1407e-11. Values are in `test_rt_budget.json`.
The outgoing energy exceeds the incident energy because the prescribed hot
electrons heat the radiation; the energy identity includes this exchange.
The run completed successfully in 26.5 wall seconds (including kernel work).

`comparison_metrics.json` records per-angle normalized spectral
L1 differences (integral of the absolute difference of unit-integral spectra),
0.27--1.06%, and limb flux ratios. These are comparison diagnostics, not proofs
of overall model accuracy or depth-resolution convergence.

An additional default eight-angle comparison is retained in
`results/eaad0101/compPS_comparison/`. At its non-native compPS angles, the limb
law differs more. This difference has not been separated into angular-grid and
compPS interpolation contributions; the ten-angle figure avoids changing the
angular nodes relative to the original paper benchmark.
