# Extraction records

`logxi1/` and `logxi2/` preserve the original per-cell Cloudy extraction
code and analysis scripts from the superseded standalone figure packages.
Their hard-coded runtime paths refer to the retained raw directories:

- `results/fe_lowion_emissivity_9e14af6b_20261010`
- `results/fe_lowion_emissivity_effac594_20261009`

These historical extraction records document n_parent, Gamma_K, yields,
transition selection, and cell-width integration. They are not required
for replotting the combined figure. Historical file-protection inventories
can reference superseded figures that have now been intentionally removed.
The current plot is reproduced with `python plot_kalpha_production.py` in
the combined package, using its included data and current source manifest.
