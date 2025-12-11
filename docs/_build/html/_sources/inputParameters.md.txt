(InputParameters)=
# Input Parameters

## Summary

Here we list the input parameters summary

[**Grids**](Grids)
| Parameter Name | Type    | Default | Description |
| :---           | :---:   | :---    | :--- |
| [`debug`](debugmodel)           | int   | 0   | Debug switch. |
| [`model`](debugmodel)           | int   | 0   | Reflection model or pure scattering. |
| [`nmaxmain`](ITERATIONstep)     | int  | 15       | Number of maximum main iteration steps. |
| [`nmaxrte`](ITERATIONstep)    | int  | 200       | Number of maximum radiative transfer iteration steps. |
| [`nfrt`](e,m,ngrids)           | int | 5000     | Number of energy bins. |
| [`nmurt`](e,m,ngrids)            | int   | 10   | Number of angle grids. |
| [`ndrt`](e,m,ngrids)            | int   | 200   | Number of depth grids. |

[**Plasma**](Plasma)
| Parameter Name | Type    | Default | Description |
| :---           | :---:   | :---    | :--- |
| [`it_temp`](temp)           | int   | 99   | Maximum iteration steps for thermal equilibrium. |
| [`temp_gas`](temp)           | float   | 1e8   | Initial gas temperature. |
| [`temp_gas_unit`](temp)          | str   | 'k'   | Unit of initial gas temperature. |
| [`nh`](nh)          | float   | 1.0e15   | Hydrogen density. |
| [`zeta`](ioniza)           | float   | 3.0  | Log10 of ionization parameter. |

[**Illumination**](illumination) 
| Parameter Name | Type    | Default | Description |
| :---           | :---:   | :---    | :--- |
| [`inci_type`](illumination)           | str   | 'powerlaw'   | Type of corona illumination. |
| [`gamma`](illumination)            | float   | 2.0   | Photon index or Blackbody temperature in eV. |
| [`hcut`](illumination)            | float   | 300e3   | High energy cut-off for cut-off power law or electrons temperture. |
| [`ktbb_nth`](illumination)            | float   | 0.1e3   | Blackbody temperature when incident type is Comptonization source. |
| [`inmu`](illumination)            | float   | 0.7  | Incidence angle. |
| [`sbot`](illumination)            | int   | 0  | Switch of bottom illumination [Thermal disk radiation]. |
| [`ktbb`](illumination)            | float   | 350.0  | Temperature of Thermal disk radiation. |
| [`fx_frac`](illumination)           | float   | 1.0  | Flux of top illumiantion over total illumination. |
| `inci_file` | str | 'incident/nthcomp_g2.0_t60.txt' | Incident spectrum file name [Necessary for type of incident is "file"]. |

[**Elemental Abundances (Relative to Hydrogen)**](abundances)
| Parameter Name | Type    | Default | Description |
| :---           | :---:   | :---    | :--- |
| `h` | float | 1.0 | Abundance of Hydrogen. This is the reference element. |
| `he` | float | 1.0 | Abundance of Helium relative to Hydrogen. |
| `li` | float | 0.0 | Abundance of Lithium relative to Hydrogen. |
| `be` | float | 0.0 | Abundance of Beryllium relative to Hydrogen. |
| `ba` | float | 0.0 | Abundance of Barium relative to Hydrogen. |
| `c` | float | 1.0 | Abundance of Carbon relative to Hydrogen. |
| `n` | float | 1.0 | Abundance of Nitrogen relative to Hydrogen. |
| `o` | float | 1.0 | Abundance of Oxygen relative to Hydrogen. |
| `f` | float | 0.0 | Abundance of Fluorine relative to Hydrogen. |
| `ne` | float | 1.0 | Abundance of Neon relative to Hydrogen. |
| `na` | float | 0.0 | Abundance of Sodium relative to Hydrogen. |
| `mg` | float | 1.0 | Abundance of Magnesium relative to Hydrogen. |
| `al` | float | 0.0 | Abundance of Aluminum relative to Hydrogen. |
| `si` | float | 1.0 | Abundance of Silicon relative to Hydrogen. |
| `p` | float | 0.0 | Abundance of Phosphorus relative to Hydrogen. |
| `s` | float | 1.0 | Abundance of Sulfur relative to Hydrogen. |
| `cl` | float | 0.0 | Abundance of Chlorine relative to Hydrogen. |
| `ar` | float | 1.0 | Abundance of Argon relative to Hydrogen. |
| `k` | float | 0.0 | Abundance of Potassium relative to Hydrogen. |
| `ca` | float | 1.0 | Abundance of Calcium relative to Hydrogen. |
| `sca` | float | 0.0 | Abundance of Scandium relative to Hydrogen. |
| `ti` | float | 0.0 | Abundance of Titanium relative to Hydrogen. |
| `va` | float | 0.0 | Abundance of Vanadium relative to Hydrogen. |
| `cr` | float | 0.0 | Abundance of Chromium relative to Hydrogen. |
| `mn` | float | 0.0 | Abundance of Manganese relative to Hydrogen. |
| `fe` | float | 1.0 | Abundance of Iron relative to Hydrogen. |
| `co` | float | 0.0 | Abundance of Cobalt relative to Hydrogen. |
| `ni` | float | 1.0 | Abundance of Nickel relative to Hydrogen. |
| `cu` | float | 0.0 | Abundance of Copper relative to Hydrogen. |
| `zn` | float | 0.0 | Abundance of Zinc relative to Hydrogen. |

**Path Group**
| Parameter Name | Type    | Default | Description |
| :---           | :---:   | :---    | :--- |
| `kernelpath` | str | 'kernel_exact5000' | Compton scattering Redistribution function file, the energy bins need to be as same in model. |
| `dataenv` | str | 'data' | Atomic database and Compton heating-cooling file. |
| `op_spec` | str | spec.dat | Output spectrum file. |
| `op_temp` | str | temp.dat | Output temperature file. |
| `op_inte` | str | inte.dat | Output intensity file. |
| `op_log` | str | model.log | Output log file. |
| `op_emis` | str | emis.dat | Output emissivity file (not used for now). |
| `op_abund` | str | abud.fits | Output abundance file. |

## Detailed Description

(Grids)=
### Grids 

(debugmodel)=
**debug,model**

DAO support the pure scattering model with constant gas temperature along the vertical
depth when `model` set at 1, and normal reflection model when `model` set at 0.

Debug mode is opened when `debug` switch set at 1. Energy, angle and depth grids and other useful inforamtion will be printed to a file when it’s open.

(ITERATIONstep)=
**nmaxmain,nmaxrte**

`nmaxrte` controls the maximum iteration steps for radiative transfer. Considering the
calculation time and convergence, we set 200 as default steps.

`nmaxmain` controls the maximum iteration steps between radiative transfer and XSTAR. We set 15 as default value


(e,m,ngrids)=
**nfrt,nmurt,nfrt**

These parameters define the resolution of depth, energy, and angle grids. We don’t
recommend user reduce the resolution of depth and energies. The photons are more
easier to be scattered in its own energy and the line emission haven’t enough resolution
when `nfrt` is small.

(plasma)=
### Plasma

(temp)=
**it_temp,temp_gas,temp_gas_unit**

XSTAR need number of iterations to get thermal equilibrium, `it-temp` is same with `nlimd`
in XSTAR, we don't recommend user change this parameter. `temp-gas`, `temp-gas-unit` define the initial gas temperature and its unit, they are only important when you run constant-temperature pure scattering model. 

(nh)=
**nh**

Hydrogen density in cm−3. The DAO version 1.0 only support for density less than 10^18
cm−3. Specify atom data table will be needed when density large than that value {cite:p}`2021ApJ...908...94K`.

(ioniza)=
**zeta**

Initial value of the log (base 10) of the model ionization parameter at the innermost shell. The {cite:p}`1969ApJ...156..943T` form is used: 

```{math}
:label: eq-xi
\xi = \frac{4\pi F_x}{n_h}
```

This value will be re-calculated after radiation field update.


(illumination)=
### Illumination

There are many parameters control the illumination, and different meaning of each parameters when we assume different type of illumination. In this section, I'll sort the parameters by differnt incident type.

#### Bottom blackbody

In `DAO`, the radiation from the bottom of the disk is modeled as a blackbody source. The parameter `sbot` controls this boundary condition: if `sbot=1`, bottom illumination is enabled; otherwise, the bottom boundary is assumed to be zero.

When `sbot=1`, `ktbb` defines the blackbody temperature, and `fx_frac` determines the partition of the total flux between the top and bottom surfaces. First, the total ionizing flux (Fx) is calculated via Eq. {eq}`eq-xi`. The fluxes for the top and bottom boundaries are then derived as:

```{math}
F_\mathrm{top} = F_x \times \mathrm{fx\_frac}, \quad F_\mathrm{bot} = F_x \times (1 - \mathrm{fx\_frac})
```

#### Power law or Cut-off power law

This formulation applies when `inci_type` is set to `powlaw` or `cutoff`. The spectral shape is defined as:

```{math}
F(E) = A \times E^{-\Gamma+1} \exp\left(-E/E_{cut}\right) , \quad \mathrm{or} \quad F(E) = A \times E^{-\Gamma+1} 
```

where: 

* `Gamma` ($\Gamma$): The photon index.

* `hcut` ($E_{cut}$): The high-energy cutoff.

* A: normalization factor, defined by Eq.{eq}`eq-xi`

Currently, the model does not implement a low-energy cutoff.

#### Blackbody

This formulation applies when `inci_type` is set to `blackbody`. The spectral shape is defined as a blackbody with temperature set at `Gamma` [eV].

#### nthcomp

This formulation applies when `inci_type` is set to `nthcomp`. The spectral shape is defined by {cite:t}`1996MNRAS.283..193Z`.

* `Gamma`: Photon index.

* `hcut`: Electron temperatures $kT_e$ [eV].

* `ktbb_nth`: Blackbody temperatures $kT_{bb}$ [eV].

#### comptt
```{versionadded} 1.1.3
This incident spectrum type was introduced in version 1.1.3. 
```
This formulation applies when `inci_type` is set to `comptt`. The spectral shape is defined by {cite:t}`1994ApJ...434..570T`. This configuration assumes a disk geometry with the optical depth fixed at $\tau = 1.0$.

* `Gamma`: Wien temperature [eV].

* `hcut`: Plasma temperature [eV].


#### file

This mode is activated when `inci_type` is set to `file`. In this configuration, the spectral shape is determined by an external user-provided file.

Consequently, the `inci_file` parameter is **mandatory** and must specify a valid path to your spectrum file.

**File Format Specification**

The spectrum file must consist of a **3-line header** followed by the data rows (total $3+n$ rows).

**1. Header Lines (First 3 rows)**

* **Row 1**: Number of energy bins ($n$).
* **Row 2**: Spectrum Unit Flag.
    * `0`: $\mathrm{erg \cdot cm^{-2} \cdot s^{-1} \cdot erg^{-1}}$
    * `1`: $\mathrm{photons \cdot cm^{-2} \cdot s^{-1} \cdot erg^{-1}}$
    * `2`: $\mathrm{erg^2 \cdot cm^{-2} \cdot s^{-1} \cdot erg^{-1}}$ ($E F_E$)
* **Row 3**: Energy Unit Flag.
    * `0`: eV
    * `1`: keV

**2. Data Block**

Following the header, the file must contain **2 columns** (Energy, Flux) for the remaining $n$ rows.

(abundances)=
### Abundances

Atomic abundances for elements H through Zn are initialized using a predefined table. The adopted values are based on {cite:t}`1996ASPC...99..117G`.

The user can change the abundance of individual elements, relative to that table.