# Figure 2 — Compton redistribution kernel slices

**Author:** Yimin Huang · Fudan University; University of Bristol · huangym23@m.fudan.edu.cn

2-D slices of the Compton redistribution kernel `K(μ_out, E_out)` at fixed
`(E_in, μ_in, T)`, for two input energies: Fe Kα (6.4 keV) and the
Compton-hump band (40 keV). Both at `T = 1e8 K`, `μ_in = 0.96`.

## Files
- `plot_kernel_redist.py` — reads the two `.dat` slices and writes the figure.
- `kernel_redist2d_FeKa.dat`, `kernel_redist2d_hump.dat` — the kernel slices.
- `kernel_redist_FeKa_hump.{png,pdf}` — the figure.
- `dump_kernel_slice.cpp` — tool that regenerates the `.dat` slices (optional).

## Run
```bash
python plot_kernel_redist.py
```
No arguments needed.

## Regenerate the data (optional)
`dump_kernel_slice.cpp` depends on the RT solver's `source/compton_kernel.cpp`
and `source/compton_cross_section.cpp`. Copy it into `source/` and build from
the repo root (`make dump_kernel_slice`), then:
```bash
./dump_kernel_slice -Ein 6.4e3 -T 1e8 -mu_in 0.96 -Eo_lo 1 -Eo_hi 500 \
                    -NMU 181 -NEO 240 -out kernel_redist2d_FeKa.dat
./dump_kernel_slice -Ein 4e4  -T 1e8 -mu_in 0.96 -Eo_lo 1 -Eo_hi 500 \
                    -NMU 181 -NEO 240 -out kernel_redist2d_hump.dat
```
