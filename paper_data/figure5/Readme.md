# Figure 5 — Iron ion-fraction vs Thomson optical depth

**Author:** Yimin Huang · Fudan University · huangym23@m.fudan.edu.cn (hyimin0924@gmail.com)

Iron ionization fraction (Fe I … Fe XXVII) vs Thomson optical depth τ_T,
shown as a colormap for two model versions:
- **(a)** DAOv1.0 (run `9bd2d8c227_C3`).
- **(b)** DAOv2.0 (run `c63d0c48`).

## Files
- `abund_9bd2d8c227_C3.fits` — DAOv1.0 iron abundances (FITS table, `fe_*` columns).
- `temp_9bd2d8c227_C3.dat` — DAOv1.0 depth/temperature table (τ from last 200 rows).
- `profile_c63d0c48_iter018.dat` — DAOv2.0 depth profile (τ_T in column 2).
- `c63d0c48<N>.iron` — DAOv2.0 iron fractions, one file per depth zone N
  (last line = latest iteration).
- `compare_iron_frac.py`, `compare_iron_frac.{png,pdf}`.

## Run
```bash
python compare_iron_frac.py
```
No arguments needed. Requires `astropy` (FITS) and `pandas`.
