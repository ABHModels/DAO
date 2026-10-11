Native Cloudy/XSTAR comparison at logxi=3, NH=1e21 cm^-2, vturbi=200 km/s.
Densities: nH=1e15 and 1e18 cm^-3. nthcomp: Gamma=1.8, kTe=60 keV, kTbb=0.01 keV.

Reproduce both figures (NumPy and Matplotlib only):
  python paper_data/figure_cloudy_xstar_nh15_nh18/plot_comparison.py
One density only: append --density 15 or --density 18.

Each lognH directory contains data/, figures/, inputs/cloudy/, inputs/xstar/, and provenance.json.
Native comparison data are sufficient to reproduce the plots; raw runtime outputs remain in results/.
All cloudy.in files, local meshes, abundance tables, xray12.ini files, and illuminating spectra are included.
XSTAR run.sh records the complete explicit parameters; parameters.json is a readable duplicate.
The Cloudy original inputs request extra diagnostic outputs; those outputs are not needed for plotting.

Panels: transmittance, reflectance, normalized reflected O VIII zoom, transmitted O VIII zoom.
Ratios are XSTAR/Cloudy on identical native bin centres. No smoothing or continuum subtraction.
Cloudy reflectance is inward gas continuum plus inward lines, excluding reflected incident continuum.
XSTAR reflectance uses emit_inward. Both are divided by their incident spectrum.
Only the reflected zoom is normalized: R(E)/mean[R(0.55 keV), R(0.75 keV)], independently per code.
Only the normalization endpoints are interpolated. The transmitted zoom is not additionally normalized.
First panel gives nH and back-face temperatures; second gives illuminated-face temperatures.

To rerun models, initialize installed HEASoft and set CLOUDY_DATA_DIR to installed Cloudy data.
Optional CLOUDY_EXECUTABLE selects Cloudy; XSTAR_EXECUTABLE selects XSTAR (default xstar).
Run scripts in fresh copies of inputs directories to avoid overwriting archived inputs.
Local continuum_mesh.ini takes precedence; no global Cloudy mesh changes are required.

The nH=1e18 result used an independently compiled XSTAR 2.59j driver with only tinf initialized
before its first use (0.099d0, matching legacy xstarsub.f). niter=99; final |H-C|/(H+C) about 7.5%.
The user accepted this residual for comparison. This result is not exact thermal equilibrium.
Installed libraries and atomic database were unchanged. nH=1e15 used installed XSTAR, niter=10.
Neither the modified XSTAR executable nor its source/build files are included here.
See XSTAR_INITIALIZATION_NOTE.txt for the exact change, source line numbers, and validation.
