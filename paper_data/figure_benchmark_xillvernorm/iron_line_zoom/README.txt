Fe K zoom of the current benchmark (log xi = 1, 2, 3)

Inputs are the curated snapshots in the parent folder:
  1: DAO 9e14af6b, final iteration 63; reference hash 1e4178a5
  2: DAO effac594, final iteration 35; reference hash 603b2ef4
  3: DAO a518ab75, final iteration 33; reference hash feb16067
No simulation or reference model was rerun. The broadband benchmark
figure and its plotting script are unchanged.

Run:
  python plot_iron_line_zoom.py
Add --show for an interactive window. Optional --anchors LOW HIGH uses
different normalization endpoints in keV and overwrites this
zoom's output, leaving the parent benchmark unchanged.
The three panels use linear flux axes in a single row.

Deliverables:
  output/pdf/iron_line_zoom_logxi123.pdf (vector figure, embedded fonts)
  output/iron_line_zoom_logxi123.png (600 dpi)
  iron_line_zoom_caption.tex (red LaTeX caption)
  data/logxi*_*.dat (the plotted spectra and ratio data)
  plot_manifest.json (normalizations, source hashes, numerical checks)

Normalization:
  Fhat(E) = F_E(E) / [(F_E(5 keV) + F_E(8 keV))/2].
Thus the arithmetic mean of the two endpoint flux densities is unity
for every model separately. The individual endpoint fluxes need not
both coincide. This replaces the parent's 20-50 keV normalization.
Reference: Walton et al. 2025, MNRAS 543, 2633, Figure 10:
  https://arxiv.org/html/2509.13411v1#S3.F10
  https://doi.org/10.1093/mnras/staf1545
Walton's oxygen comparison used 0.55 and 0.75 keV; 5 and 8 keV here
are our chosen Fe-band endpoints, not energies specified by Walton.

Flux conventions match the existing benchmark: DAO is the upper-face
outward normal flux 2*pi*sum(w_i*mu_i*I_i), mu_i>0. The energy axes are
converted to keV. Reference E*F_E columns are divided by energy to
recover F_E before endpoint normalization. Their overall input units
cancel in this normalization.

Complete-spectrum comparison:
The main panels and inset ratios use the complete endpoint-normalized
spectra, including the continuum. No continuum fitting or subtraction,
line-power normalization or equivalent-width measurement is performed.
The saved model data contain two columns: E_keV and Fhat_E.

Resolution and ratios:
Main panels connect native data samples linearly; no Gaussian smoothing
or relativistic/instrumental broadening is added. Native peak heights
and widths therefore reflect different model energy resolutions.
Inset ratios use the complete normalized flux averaged
into bins with geometric-midpoint edges around the reflionx samples
between 6.15 and 7.10 keV. Integrate the piecewise-linear input spectra
over each bin without continuum subtraction.
This avoids comparing peaks evaluated at different sample locations
and preserves the integrated input flux over the selected bins; it
does not remove all differences in intrinsic spectral resolution.
The ratio is Fhat_reference/Fhat_DAO, with DAO/DAO=1. Bins with
non-finite or non-positive averaged DAO flux are excluded. Validity is
recorded explicitly in the saved ratio data. These are spectral-shape
ratios under the stated normalization, not absolute-flux ratios.

The figure uses one row and three columns. The inset positions vary
to avoid the dominant emission peak. The PDF is intended for full
page width in a two-column manuscript.
