# Local Fe K-alpha production at log xi = 2

This is the first proposed diagnostic only: calculate local line production
and integrate it over the slab. No temperature search or radiative-transfer
solve is performed, and no boundary illumination is changed.

## Accepted atmosphere

Source: `results/effac594`, converged outer iteration 35. Reproduce the
constant-temperature Cloudy pass using temperatures from
`thermal_budget_iter35.dat` and J from `moments_iter034.dat`. This is the
preceding-outer J actually used to produce the saved accepted atomic state;
using iteration 35 J would instead perform a new atomic update. Hydrogen
density is 1e15 cm^-3, iron scaling is 1, with the original 48 depth cells.

The diagnostic includes the unchanged production atomic interface and
appends exports after extraction. It compiles current source files into a
separate executable. All changes and outputs stay in this directory.
The replayed iron fractions agree with the saved fractions to a maximum
absolute difference of 6.33e-8. The 363 protected original files, including
production source, the source run and prior figures, were unchanged.

## What the numbers mean

For each inner-shell fluorescence entry, export the absorbing parent ion
density n, K-shell photoionization rate Gamma, photon yield omega and energy.
The all-direction local power is epsilon = n Gamma omega E.

Select Fe K-shell vacancies and line types 1/2 (K-alpha doublet) in the
actual Cloudy `mewe_fluor.dat`; exclude types 3/4 (K-beta). Also collect the
bound-bound Fe K-alpha complex between 6.2 and 7.0 keV using actual
transition energies, excluding bookkeeping lines and blends as in production.
The denominator includes fluorescence plus these bound-bound lines before
DAO survival correction. Post-survival material-source totals are supplied
separately; neither quantity is the emergent spectrum.

Fluorescence is grouped by the ion BEFORE K-shell photoionization, which is
the population in the ion-fraction plot. The emitting ion can have a higher
charge. Both parent and emitting stages are exported. Bound-bound entries
are instead grouped by their emitting ion. This distinction must be kept
when discussing which ion populations are absent from REFLIONX.

Integrate as sum(epsilon_d * dz_d), using the actual nonuniform cell widths.
The result is generated power per slab area in erg cm^-2 s^-1, integrated
over all emission directions. It is NOT an upper-face flux or an equivalent
width. Local epsilon has units erg cm^-3 s^-1.

Imported Rydberg energies use DAO's existing 13.6058 eV/Ryd conversion,
which differs from Cloudy's database conversion by about 8 ppm. Both line
energy and assigned RT-bin energy are retained, and totals using the
production RT-bin energies are also available in `analysis.json`.

## Results

- All Fe K-alpha generated power per area: 8.7441e13 erg cm^-2 s^-1.
- Fe I-V parent fluorescence: 3.0013e12 erg cm^-2 s^-1, or 3.4324%.
- Fe VI-XVI parent fluorescence: 85.6261% of all Fe K-alpha production.
- Fe XVII-XXII parent fluorescence: 10.9415%.
- Bound-bound Fe K-alpha production: 3.4483e7 erg cm^-2 s^-1.
- Fe I/II contributions are zero; the Fe I-V contribution is almost all Fe V.
- 92.23% of Fe I-V K-alpha production occurs in cells at tau_ref > 1.
- The maximum local Fe I-V contribution is 46.16%, but it is confined to
  layers whose absolute K-alpha emissivity is much lower than at the surface.

These are production fractions, not escaping fractions. They do not quantify
the change in the observed iron line or establish a unique explanation of the
DAO-REFLIONX difference. No source-removal or zero-boundary experiment is part
of this diagnostic.

## Files

- `depth_power.dat`: local powers, depth-weighted contributions and fractions.
- `ion_power.dat`: slab-integrated powers by ion stage.
- `fluorescence_records.dat`: individual source factors, identities and rates.
- `cells/`: individual per-depth fluorescence and bound-bound exports.
- `bound_bound_survival.dat`: original DAO survival factors for source accounting.
- `local_Kalpha_production.png`: local production and Fe I-V power fraction.
- `analysis.json`, `provenance.json`: numerical summary and input hashes.
- `protected_files_before.json`: hashes of original files preserved by this task.

## Reproduce

From the DAO repository root:

```sh
make -f Makefile -f results/fe_lowion_emissivity_effac594_20261009/diagnostic.mk results/fe_lowion_emissivity_effac594_20261009/diagnostic
/usr/bin/env DYLD_LIBRARY_PATH=/Users/xe26734/Downloads/heasoft-6.37.1/aarch64-apple-darwin25.6.0/lib:/opt/homebrew/lib DAO_CLOUDY_WORKERS=4 results/fe_lowion_emissivity_effac594_20261009/diagnostic
/Applications/anaconda3/bin/python results/fe_lowion_emissivity_effac594_20261009/analyze.py
```
