# Figure 3 — DAO vs `reflionx` vs `xillvercp` over a log ξ scan

**Author:** Yimin Huang
**Affiliation:** Fudan University
**Email:** huangym23@m.fudan.edu.cn (hyimin0924@gmail.com)

Data and script to reproduce Figure 3: emergent reflection spectra from
DAOv2.0 compared against `reflionx` and `xillvercp` at log ξ = 1, 2, 3.

## Files

| Hash       | log ξ | Panel |
|------------|------:|-------|
| `c63d0c48` | 1     | a     |
| `8b877447` | 2     | b     |
| `50bd6dd6` | 3     | c     |

For each hash:

- `dao_<hash>.dat` — DAO emergent intensity (E, incident, I at 8 μ nodes)
- `spectra_<hash>.dat` — `reflionx` spectrum (E [keV], E·F_E)
- `xillver_<hash>.dat` — `xillvercp` spectrum (E [keV], E·F_E)
- `pa<hash>.json` — run parameters

All three runs share: `nthcomp` corona, Γ = 2, kT_e = 60 keV,
kT_bb = 0.01 keV, kT_disk = 0.35 keV, log n_H = 15, μ_inc = cos 45°.

`emergent_iter020.dat` is a copy of the log ξ = 3 DAO output.

## Plot

```bash
python plot_compare_reflionx.py
```

Note: the script has absolute paths hardcoded at the top — edit them to
point at the files in this folder before running.
