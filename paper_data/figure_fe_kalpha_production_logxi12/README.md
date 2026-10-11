# Local Fe K-alpha production: log xi = 1 and 2

Combined publication figure from the existing accepted atmospheres:

- log xi = 1: `9e14af6b`, iteration 63; solid curves.
- log xi = 2: `effac594`, iteration 35; dashed curves.

Dark gray denotes all Fe K-alpha production, purple Fe I-V, and teal Fe XVII-XXII.
The upper panel shows absolute local emissivity in erg cm^-3 s^-1, without
rescaling either model. The lower panel divides each group's emissivity by
the local all-Fe K-alpha emissivity of the SAME log xi model, in percent.
The upper panel omits exact zeros only because its ordinate is logarithmic;
the saved numbers and lower panel retain them. No smoothing or new model
calculation is applied.

Both depth axes use the reference Thomson optical depth, defined in the
source atmospheres with n_e,ref = 1.21 n_H. The figure retains the previous
single-column size, no title, and legends inside the panels.

## Reproduce

```sh
python plot_kalpha_production.py
```

Requires only NumPy and Matplotlib; add `--show` for an interactive window.
The package contains its own copied input tables and verifies their hashes
before plotting. It also verifies the width-weighted integrated powers
against the original extraction summaries. Outputs:

- `output/pdf/local_Kalpha_production_logxi12.pdf`: vector figure.
- `output/local_Kalpha_production_logxi12.png`: 600 dpi preview.
- `caption.tex`: red LaTeX caption.
- `plot_values.json`: definitions and integrated powers for both models.
- `source_snapshot.json`: source locations and SHA-256 hashes.
- `artifact_hashes.json`: hashes of package artifacts.

## Production definition

Fluorescence is grouped by the parent ion BEFORE K-shell photoionization.
The source is n_parent * Gamma_K * photon_yield * photon_energy, with K-shell
vacancies and database line types 1/2 (K-alpha), excluding K-beta. All-Fe
production also includes bound-bound Fe lines at 6.2-7.0 keV before DAO's
line-survival correction. The detailed extraction code and diagnostics are
consolidated in `extraction/logxi1/`, `extraction/logxi2/`, and the
corresponding `data/` subdirectories. See `extraction/README.md` for the
original runtime locations.

Whole-slab fractions use sum(epsilon_group * dz) / sum(epsilon_all * dz),
with the actual nonuniform cell widths. Fe I-V supplies 89.400462% at log
xi = 1 and 3.432415% at log xi = 2. Fe XVII-XXII supplies 0.0000327290%
and 10.941452%, respectively. These are generated powers, not emergent fluxes.
The standalone figure packages have been removed; this combined version
is the retained publication figure. Original model and extraction runs
remain in `results/`.
