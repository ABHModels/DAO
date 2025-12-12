(outputfile)=
# Output Files

## Summary

Here we list the output files summary. We use the example in [Sample Output](sampleoutput)
| File Name | description |
| :---      | :-------- |
|`[out/User_Guide/spec_b912ee2e38_User_Guide.dat](spec.dat)`|energy, top illumination, bottom illumination and angle averaged intensity at top and bottom layer|
|`[out/User_Guide/temp_b912ee2e38_User_Guide.dat](temp.dat)`|Thomson optical depth, temperature|
|[`out/User_Guide/inte_b912ee2e38_User_Guide.dat`](inte.dat)|energy, intensity for 10 angle bins|
|[`out/User_Guide/model_b912ee2e38_User_Guide.log`](log.log)|model log|
|[`out/User_Guide/abund_b912ee2e38_User_Guide.fits`](abund.fits)|same with XSTAR [xout_abund1.fits](https://heasarc.gsfc.nasa.gov/docs/software/xstar/docs/sphinx/xstardoc/docs/build/html/xstaroutput.html#the-abundances-data-file-xout-abund1-fits)|


## Detailed Description

(spec.dat)=
### spec.dat
* **Column 1**: Energy grid values [ev].
* **Column 2**: The external illuminaiton at upper surface [erg/cm2/s/erg/ster].
* **Column 3**: The external illuminaiton at lower surface [erg/cm2/s/erg/ster].
* **Column 4**: The emergent angle-averaged spectrum at lower surface [erg/cm2/s/erg].
* **Column 5**: The emergent angle-averaged spectrum at upper surface [erg/cm2/s/erg].

We record these data for each iterations.

(temp.dat)=
### temp.dat
* **Column 1**: Thomson optical depth.
* **Column 2**: Plasma temperature at each depth [K].

We record these data for each iterations.

(inte.dat)=
### inte.dat
* **Column 1**: Energy grid values [ev].
* **Column 2 through 111**: The emergent specific intensity [erg/cm2/s/erg/ster] for each of 10 angle bins [0.05 to 0.95 with linear grids]. Column 2 contains the intensity for the first angle ($\mu_1$), Column 3 for the second angle ($\mu_2$), and so on. 

We record these data for each iterations.

(log.log)=
### model.log
Model log file

(abund.fits)=
### abund.fits
This ASCII FITS file contains ion abundances as well as heating and cooling rates. Only ions with fractional abundance (relative to its parent element) greater than $10^{-10}$ (relative to its parent element) are recorded. The elements are ordered by increasing nuclear charge, ions by increasing free charge. See [The Abundances Data File: xout_abund1.fits](https://heasarc.gsfc.nasa.gov/docs/software/xstar/docs/sphinx/xstardoc/docs/build/html/xstaroutput.html#the-abundances-data-file-xout-abund1-fits) for more details.

We record these data for final results.