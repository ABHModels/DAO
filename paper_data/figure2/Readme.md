# Figure 2 — Compton redistribution kernel slices

**Author:** Yimin Huang
**Affiliation:** Fudan University
**Email:** huangym23@m.fudan.edu.cn (hyimin0924@gmail.com)

2-D slices of the Compton redistribution kernel `K(μ_out, E_out)`
at fixed `(E_in, μ_in, T)`, for two input energies:
Fe Kα (6.4 keV) and the Compton-hump source band (40 keV).

## Files

- `dump_kernel_slice.cpp` — tool that writes a kernel slice.
- `kernel_redist2d_FeKa.dat` — slice at `E_in = 6.4 keV`, `T = 1e8 K`, `μ_in = 0.96`.
- `kernel_redist2d_hump.dat` — slice at `E_in = 40 keV`, same `T`, `μ_in`.
- `plot_kernel_redist.py` — produces the figure from the two `.dat` files.

## Build

`dump_kernel_slice.cpp` depends on `source/compton_kernel.cpp` and
`source/compton_cross_section.cpp`. Copy it into `source/`, then
build from the repo root:

```bash
cp paper_data/figure2/dump_kernel_slice.cpp source/
make dump_kernel_slice
```

## Regenerate the data

```bash
./dump_kernel_slice -Ein 6.4e3 -T 1e8 -mu_in 0.96 \
                    -Eo_lo 1 -Eo_hi 500 -NMU 181 -NEO 240 \
                    -out kernel_redist2d_FeKa.dat

./dump_kernel_slice -Ein 4e4  -T 1e8 -mu_in 0.96 \
                    -Eo_lo 1 -Eo_hi 500 -NMU 181 -NEO 240 \
                    -out kernel_redist2d_hump.dat
```

## Plot

```bash
python plot_kernel_redist.py
```