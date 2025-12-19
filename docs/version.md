(versionpage)=
# Version

## v1.1

Standard model framework

## v1.1.2

* Add thermally comptonized continuum `nthcomp` as a type of incident spectrum. When incident type is `nthcomp`, `hcut` is used to represent corona temperature in eV. The seed photons' temperature is fixed at 0.1 keV.

* Add new atom database for hydrogen density less than 10$^{22}$. To use high density, you must denote `dataenv` as a right folder which include atom database for $n_h\le10^{18}$ (*dataenv/len18*), $10^{18}<n_e\le10^{19}$ (*dataenv/n19*), $10^{19}<n_e\le10^{20}$ (*dataenv/n20*), $10^{20}<n_e\le10^{21}$ (*dataenv/n21*),  $10^{21}<n_e\le10^{22}$ (*dataenv/n22*).

## v1.1.3

* Add `comptt` as a type of incident spectrum
* When incident spectrum is `nthcomp`, the blackbody temperature is a free parameter now.
* Fix the normalization error for bottom disk radiation. **[Contributor: Lin,Gao, Fudan University]**

## v1.1.4

* Add low energy cutoff for `powerlaw` and `cutoff` incident. The low cutoff energy is fixed at 0.1keV
$$
    F(E)=E^{-\Gamma+1}\times\exp{-\frac{E}{E_{hcut}}}\times\exp{-\frac{E_{lcut}}{E}},\quad E_{lcut}=0.1 keV
$$