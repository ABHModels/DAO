# Fe XXV survival correction at log xi = 3

Three-curve comparison for Fig. `fe_survival_comparison`: XILLVER, original
DAO, and DAO with only the Fe XXV w and x survival factors set to one.
Source atmosphere: `a518ab75`, accepted outer iteration 33.

## Reproduce

```sh
python plot_survival_comparison.py
```

The plotting script uses only the curated files in `data/` and needs NumPy
and Matplotlib. It does not run Cloudy, a temperature search or radiative
transfer. It can be invoked from any directory; paths resolve relative to
the script.

- `output/pdf/fe_xxv_model_comparison.pdf`: vector paper figure.
- `output/fe_xxv_model_comparison.png`: 350 dpi preview.
- `output/fe_xxv_model_comparison_manifest.json`: plotting conventions and hashes.
- `caption.tex`: red LaTeX caption with the paper figure label.
- `data/`: full upper-face spectra,
  benchmark reference spectra, source data, validation and normalization.
- `source_snapshot.json`: original file locations, checksums and unit conversions.
- `experiment/`: C++ driver and build instructions for the two complete spectra.

## Meaning of the control

The strongest DAO peak contains the Fe XXV resonance line w at 6.700 keV
and magnetic-quadrupole line x at 6.682 keV. Their local surviving powers
epsilon_w P_w and epsilon_x P_x are replaced by epsilon_w and epsilon_x.
The original illumination, gas temperature, ion populations, opacities and
all other material emission stay fixed. The shared production transfer
solver iterates the Compton source. This is a fixed-atmosphere response,
not a newly converged thermal solution or an all-lines-P=1 experiment.

The figure retains the original native grids and three curves.
EW annotations have been removed; plotting does not calculate or display EWs.
Both DAO curves use the original DAO mean F_E at 5 and 8 keV as a common
normalization divisor. The XILLVER-CP reference uses its own endpoint divisor.
No extra broadening or continuum subtraction is applied in this figure.

The full spectral peak increases by 4.9852% in the fixed-atmosphere control.
The experiment driver retains its original runtime/build paths; repeating
that calculation requires the original atmosphere and dependencies. Only the
portable plotting script is needed to reproduce the figure.
