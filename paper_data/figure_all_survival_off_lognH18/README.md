# Fully re-equilibrated all-line P=1 comparison

The baseline is DAO run 5b7fe321. The comparison is run fe48cbde with the
additional DAO survival correction disabled for ALL bound-bound lines.
Temperature, ion populations, opacities, emissivities and radiation are solved
again from normal initialization with unchanged convergence criteria.
This is distinct from the fixed-atmosphere O VIII-only test.

The incident parameters are log xi=3, log nH=18, Gamma=1.8, kTe=60 keV,
kTbb=0.01 keV, reflionx flux normalization, solar Fe, 48 depth cells,
reference Thomson depth 5, angle-averaged Compton scattering.

Run `python plot_comparison.py` here to recreate the figure using only the
curated inputs (NumPy and Matplotlib required). All panels show EF_E and retain
native energy sampling without smoothing. BOTH DAO curves share the STANDARD
DAO normalization in each panel. XillverCp and reflionx retain their individual
normalizations from the previous comparison. Absolute outgoing fluxes are
stored in input/DAO_standard.dat and input/DAO_all_P1.dat, before normalization.
The units and angular-flux definition are recorded in their headers.

The line probabilities in this artificial P=1 experiment are imposed rather
than inferred from its optical depths. Continuum absorption and Compton
scattering remain active. Collisional rates in Cloudy's statistical equilibrium
are retained. Fluorescence and the existing H I Lyalpha thin-reference
normalization remain unchanged.

measurements.json records the final convergence criteria, total runtime,
absolute O VIII peak change and input hashes. Changes in the full line peak
include changes in the atmosphere and other emission; they are not an isolated
measurement of O VIII photon survival at fixed structure.

The experiment directory records the exact source patch, isolated build and
launch scripts; validation contains source/input provenance and the two-path
line-source verification. Repeating the thermal solve requires the original
DAO/Cloudy/HEASoft installation, kernel and source versions. Production source
files and the original accepted results were verified unchanged.
