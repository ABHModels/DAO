(versionpage)=
# Version

## v1.0

Standard model

## v1.1.2

* Add thermally comptonized continuum `nthcomp` as a type of incident spectrum. When incident type is `nthcomp`, `hcut` is used to represent corona temperature in eV. The seed photons' temperature is fixed at 0.1 keV.

* Add new atom database for hydrogen density less than 10$^{22}$. To use high density, you must denote \texttt{dataenv} as a right folder which include atom database for $n_h\le10^{18}$ (\texttt{dataenv/len18}), $10^{18}<n_e\le10^{19}$ (\texttt{dataenv/n19}), $10^{19}<n_e\le10^{20}$ (\texttt{dataenv/n20}), $10^{20}<n_e\le10^{21}$ (\texttt{dataenv/n21}),  $10^{21}<n_e\le10^{22}$ (\texttt{dataenv/n22}).

## v1.1.3

* Add `comptt` as a type of incident spectrum

* Fix the normalization error for bottom disk radiation. [Contributor: Lin,Gao]