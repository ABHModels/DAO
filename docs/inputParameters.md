(InputParameters)=
# Input Parameters

## Summary

Here we list the input parameters summary

**Grids**
| Parameter Name | Type    | Default | Description |
| :---           | :---:   | :---    | :--- |
| `debug`           | int   | 0   | Debug switch. |
| `nmaxmain`     | int  | 15       | Number of maximum main iteration steps. |
| `nmaxrte`    | int  | 200       | Number of maximum radiative transfer iteration steps. |
| `nfrt`           | int | 5000     | Number of energy bins. |
| `nmurt`           | int   | 10   | Number of angle grids. |
| `ndrt`           | int   | 200   | Number of depth grids. |

**Plasma**
| Parameter Name | Type    | Default | Description |
| :---           | :---:   | :---    | :--- |
| `it_temp`           | int   | 99   | Maximum iteration steps for thermal equilibrium. |
| `temp_gas`           | float   | 1e8   | Initial gas temperature. |
| `temp_gas_unit`           | str   | 'k'   | Unit of initial gas temperature. |
| `nh`           | float   | 1.0e15   | Hydrogen density. |
| `zeta`           | float   | 3.0  | Log10 of ionization parameter. |

**Illumination**
| Parameter Name | Type    | Default | Description |
| :---           | :---:   | :---    | :--- |
| `inci_type`           | str   | 'powerlaw'   | Type of corona illumination. |
| `gamma`           | float   | 2.0   | Photon index or Blackbody temperature in eV. |
| `hcut`           | float   | 300e3   | High energy cut-off for cut-off power law or electrons temperture. |
| `ktbb_nth`           | float   | 0.1e3   | Blackbody temperature when incident type is Comptonization source. |
| `inmu`           | float   | 0.7  | Incidence angle. |
| `sbot`           | int   | 0  | Switch of bottom illumination [Thermal disk radiation]. |
| `ktbb`           | float   | 350.0  | Temperature of Thermal disk radiation. |
| `fx_frac`           | float   | 1.0  | Flux of top illumiantion over total illumination. |

**Elemental Abundances (Relative to Hydrogen)**
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