# All-line survival-off thermal experiment

Baseline: converged `results/5b7fe321`, log nH=18, Gamma=1.8, log xi=3,
kTe=60 keV, kTbb=0.01 keV, reflionx incident normalization, 48 depth cells,
reference Thomson depth 5, angle-mean Compton scattering, four Cloudy workers.

New run: `results/fe48cbde`. Every bound-bound line source is epsilon/(4 pi dE),
with P=1, in both the local temperature search and full-column transfer.
The same production temperature and transfer solvers, iteration limits,
energy checks and five convergence criteria are retained. This calculation
starts from the normal boundary-field initialization and recomputes temperature,
ionization, emissivity, opacity and radiation together. It is not the earlier
fixed-atmosphere O VIII-only source control.

Continuum absorption, Compton scattering, the incident corona, atomic statistical
equilibrium and collisional rates are retained. Fluorescence is unchanged
(it already bypasses the DAO bound-bound survival factor). The existing H I
Lyalpha thin-reference correction is also retained. Disabling the additional
DAO survival factor does not disable atomic collisional de-excitation inside
Cloudy's statistical equilibrium.

The isolated `line_probabilities` routine returns P=1 for every transition.
Its beta=1, Pelec=0 and Pdest=0 are imposed experimental values; the saved
probability diagnostics must not be interpreted as predictions of physical
line escape from the new atmosphere. Trapping-destruction diagnostics are zero.

Only copies in this directory are edited. `experiment.patch` records every
change relative to the production sources. A distinct physics hash and a
refusal to overwrite an existing result directory keep the baseline separate.
The standalone executable's verification checks both line-source paths with
opaque, collisionally quenched synthetic lines, multiple lines per bin and
a zero-emission absorber, while checking that continuum, fluorescence and
opacities are preserved.

Build from the repository root (requires the original Cloudy/HEASoft installation):

```bash
/Applications/anaconda3/bin/python results/all_survival_off_5b7fe321/prepare.py
make -f Makefile -f results/all_survival_off_5b7fe321/experiment.mk results/all_survival_off_5b7fe321/maindaocl_all_P1
/Applications/anaconda3/bin/python -u results/all_survival_off_5b7fe321/run.py
```

The preparation and launch scripts deliberately refuse to overwrite the run.
For an exact record of this execution, see `launch.json`, `console.log` and
`results/fe48cbde/thermal_status.json`. A saved intermediate spectrum is not a
converged result. `collect.py` accepts only a completed run passing all five
conditions for the final three iterations, verifies all emitting lines have
P=1, and creates the curated `paper_data/figure_all_survival_off_lognH18` package.
