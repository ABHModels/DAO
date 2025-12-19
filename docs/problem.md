(problem)=
# Known Issues and Limitations

This section documents known issues and current limitations of the `DAO` model. It will be updated continuously as development progresses. Users are encouraged to report any potential discrepancies or bugs found during usage.

* **Bottom Illumination:** The current release has not yet been rigorously tested for the bottom illumination scenario. A detailed analysis and validation of this feature will be presented in a forthcoming publication.
* **High-Density Regime:** Discrepancies have been observed in the high-density version of the model, specifically at a density of $n_{\rm H}=10^{21}$ (likely $\rm cm^{-3}$). As illustrated in the figure below, artifacts appear at this density level. Preliminary analysis suggests that these errors may stem from inaccuracies in the underlying atomic data tables.

```{figure} images/highdensity_with_1e21.svg
:width: 100%
:align: center
:name: fig-highdensity_with_1e21

Spectra generated with `nthcomp` incident radiation, including the density point $n_{\rm H}=10^{21}$.
```
```{figure} images/highdensity_without_1e21.svg
:width: 100%
:align: center
:name: fig-highdensity_with_1e21

Same as above but exculde the $10^{21}$
```

