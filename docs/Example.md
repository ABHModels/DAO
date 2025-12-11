(example)=

# Example

## Effect of ionization parameters

The Python interface is particularly powerful for running a grid of models on a cluster. The example `run.sh` script, shown in Listing follows, is configured to demonstrate this by looping over 16 different values of the ionization parameter, `zeta`.

<details>
<summary>Click to expand run.sh script</summary>

```bash
PYTHON_SCRIPT="modelrun.py" # <-- Make sure this is your main Python script name
DIY="dxi5000"

# --- Execution & Model Flags ---
exec=1      # (0: Generate files only, 1: Submit with sbatch, 2: Run locally)
debug=0     # Debug switch
model=0     # (0: Reflection, 1: Pure scattering)
ity=2       # Type of angle bins and integral method

# --- Computational Grids ---
nmaxmain=15        # Number of maximum main iteration steps
nmaxrte=100        # Number of maximum radiative transfer iteration steps
ndrt=200           # Number of depth grids (tau)
nfrt=5000          # Number of energy grids [Must match kernel file]
nmurt=10           # Number of angle grids (mu)
ppemin=0.1         # Minimum energy of the grid [eV]
ppemax=9e5         # Before maximum energy of the grid [eV]
ppemax2=1e6        # Maximum energy of the grid [eV]
taumin=1e-4        # Minimum Thomson optical depth of the grid
taumaxrt=10.0      # Maximum Thomson optical depth of the grid
mumin=0.05         # Minimum mu value for the angular grid (> 0)
mumax=0.95         # Maximum mu value for the angular grid (< 1)

# --- Temperature & Equilibrium ---
it_temp=99         # Maximum iteration steps for thermal equilibrium
temp_gas=1e8       # Initial gas temperature
temp_gas_unit='k'  # Unit of the initial temperature ('k' or 'ev')

# --- Incident Spectrum Parameters ---
inci_type="file"               # Type of corona illumination ('powerlaw', 'file', 'cutoff')
gamma=1e3                           # Photon index or Blackbody temperature in eV
hcut=300e3                       # High energy cut-off for cut-off power law illumination [eV]
inmu=0.7                         # Cosine of the incident illumination angle
inci_file="incident/nthcomp_g2.4_t60.txt" # Path to incident spectrum file

# --- Bottom Illumination (Disk) ---
sbot=0              # Switch for bottom illumination (0: off, 1: on)
ktbb=350.0          # Temperature of the bottom blackbody source [eV]
fx_frac=1e-10       # Flux ratio of top illumination to total illumination

# --- Atmosphere Properties ---
zeta=3.0            # Log10 of ionization parameter (xi)
nh=1e15             # Hydrogen number density [cm^-3]
fe=1.0              # Iron abundance relative to Hydrogen

# --- Paths ---
kernelpath="kernel_exact5000" # Path to the Compton scattering redistribution file

echo "----------------------------------------"
echo "Starting model run with DIY tag: ${DIY}"
echo "----------------------------------------"

# Example loop: You can loop over any parameter.
# Here, we loop over the incident angle 'inmu'.
for zeta in 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.2 3.4 3.6 3.8 4.0
for inci_file in "incident/nthcomp_g2.4_t60.txt" "incident/nthcomp_g1.4_t60.txt"
do
    echo "Running with nh = ${zeta}"
    
    python ${PYTHON_SCRIPT} \
        --exec ${exec} \
        --debug ${debug} \
        --model ${model} \
        --ity ${ity} \
        --nmaxmain ${nmaxmain} \
        --nmaxrte ${nmaxrte} \
        --ndrt ${ndrt} \
        --nfrt ${nfrt} \
        --nmurt ${nmurt} \
        --ppemin ${ppemin} \
        --ppemax ${ppemax} \
        --ppemax2 ${ppemax2} \
        --taumin ${taumin} \
        --taumaxrt ${taumaxrt} \
        --mumin ${mumin} \
        --mumax ${mumax} \
        --it-temp ${it_temp} \
        --temp-gas ${temp_gas} \
        --temp-gas-unit "${temp_gas_unit}" \
        --inci-type "${inci_type}" \
        --gamma ${gamma} \
        --hcut ${hcut} \
        --inmu ${inmu} \
        --inci-file "${inci_file}" \
        --sbot ${sbot} \
        --ktbb ${ktbb} \
        --fx-frac ${fx_frac} \
        --zeta ${zeta} \
        --nh ${nh} \
        --fe ${fe} \
        --kernelpath "${kernelpath}" \
        --diy "${DIY}"

    echo "----------------------------------------"
done
done
```
</details>

<details>
<summary>Click to expand plot script</summary>

```python
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl 
import json
import os 
from rdata import read_spec,read_spec5000 
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from datetime import datetime, time
def normtopflux(zeta,nh,top,ener):
    fx = 10**zeta*nh/4/np.pi
    integral = np.trapz(top,ener)
    normfactor = integral/fx
    return top*normfactor

plt.style.use('seaborn-v0_8-ticks')
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = r'\usepackage{mathptmx}' 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.size'] = 20

ERGSEV = 1.60217662e-12 

# All metadata data
METADATA_DIR = 'metadata'
if not os.path.isdir(METADATA_DIR):
    raise FileNotFoundError(f"no '{METADATA_DIR}' file")
datas = os.listdir(METADATA_DIR)

# get all data
plot_data_list = []
for data_filename in datas:
    with open(os.path.join(METADATA_DIR, data_filename), 'r') as f:
        try:
            data = json.load(f)
            params = data['parameters']
            if params.get('diy_tag') == 'dxi5000':
                try:
                    ener,top,bot,inte,flux = read_spec5000(params['op_spec'])
                    task_info = {
                        'zeta' : float(params['zeta']),
                        'ener': ener, 
                        'eemo': flux,
                        'top': normtopflux(1.0, 1e15, top, ener),
                        'Afe': params['fe'],
                        'incident': params['inci_file']
                    }
                    plot_data_list.append(task_info)  
                except Exception as e:
                    continue   
        except (json.JSONDecodeError, KeyError, TypeError) as e:
            print(f"Warning: skip {data_filename}: {e}")
            continue

if not plot_data_list:
    print("no 'diy_tag: dxi' ")
    exit()

plot_data_list.sort(key=lambda item: item['zeta'])
zeta_values = [item['zeta'] for item in plot_data_list]
inci_files = ['incident/nthcomp_g1.4_t60.txt', 'incident/nthcomp_g2.4_t60.txt']

cmap = mpl.cm.get_cmap('viridis')
norm = Normalize(vmin=min(zeta_values), vmax=max(zeta_values))
sm = ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([]) 

fig, ax = plt.subplots(
    len(inci_files), 2,           
    sharex='col',
    sharey='row',
    constrained_layout=True,
    gridspec_kw={'width_ratios': [3, 1]}  
)

main_ax = ax[:, 0]
zoom_ax = ax[:, 1]

fig.supylabel(r'$E J_E \quad [\mathrm{erg}^{2}  \cdot \mathrm{cm}^{-2}  \cdot \mathrm{s}^{-1} \cdot  \mathrm{erg}^{-1}]$')

top_plotted = [False, False]
for item in plot_data_list:
    color = cmap(norm(item['zeta']))
    energ = item['ener']
    flux = item['eemo']
    top = item['top']
    
    lw = 1.0 

    if item['incident'] == inci_files[0]:
        main_ax[0].loglog(energ, flux, color=color, lw=lw)
        if not top_plotted[0]:
            axs = main_ax[0]
            axs.text(0.05, 0.95, rf'$\Gamma$ = 1.4', transform=axs.transAxes,
                     va='top', ha='left')
            main_ax[0].loglog(energ, top, color='k', lw=2,
                              linestyle='--', label='Incident')
            top_plotted[0] = True

    if item['incident'] == inci_files[1]:
        main_ax[1].loglog(energ, flux, color=color, lw=lw)
        if not top_plotted[1]:
            axs = main_ax[1]
            axs.text(0.05, 0.95, rf'$\Gamma$ = 2.4', transform=axs.transAxes,
                     va='top', ha='left')
            main_ax[1].loglog(energ, top, color='k', lw=2,
                              linestyle='--', label='Incident')
            top_plotted[1] = True

for i, axis in enumerate(main_ax):
    axis.set_xlim(0.001, 1e3)
    axis.set_ylim(10**6.5, 1e15)
    if i == 0:
        handles, labels = axis.get_legend_handles_labels()
        if handles:  
            axis.legend(loc='lower right')

main_ax[-1].set_xlabel('Energy [keV]')
xticks = [1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]
xticklabels = [r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$',
               r'$1$', r'$10^{1}$', r'$10^{2}$', r'$10^{3}$']
main_ax[-1].set_xticks(xticks)
main_ax[-1].set_xticklabels(xticklabels)


for idx in [0, 1]:
    axz = zoom_ax[idx]
    
    zoom_flux_list = []
    zoom_ener_list = []
    for item in plot_data_list:
        if item['incident'] == inci_files[idx]:
            color = cmap(norm(item['zeta']))
            ener = np.array(item['ener'])
            flux = np.array(item['eemo'])
            
            mask = (ener >= 1.0) & (ener <= 10.0)
            if mask.any():
                axz.loglog(ener[mask], flux[mask], lw=0.8, color=color)
                zoom_flux_list.append(flux[mask])
                zoom_ener_list.append(ener[mask])

    axz.set_xlim(1, 10)

    if zoom_flux_list:
        ymin = min(np.min(f) for f in zoom_flux_list)
        ymax = max(np.max(f) for f in zoom_flux_list)
        axz.set_ylim(ymin * 0.8, ymax * 1.2)


    axz.tick_params(axis='y', which='both', labelleft=False)

    axz.text(0.05, 0.95, '1–10 keV', transform=axz.transAxes,
             va='top', ha='left')


zoom_ax[-1].set_xlabel('Energy [keV]')
zoom_ax[-1].set_xticks([1e0, 1e1])
zoom_ax[-1].set_xticklabels([r'$1$', r'$10$'])
zoom_ax[-1].tick_params(axis='x', which='minor', labelbottom=False)
cbar = fig.colorbar(sm, ax=ax.ravel().tolist(), location='right', shrink=0.8)
cbar.ax.set_xlabel(r'log$\xi$')
cbar.ax.tick_params(labelsize=15)

plt.show()
```
</details>


```{figure} images/dxi.svg
:width: 100%
:align: center
:name: fig-dxi

Angle-averaged reflection spectra for various ionization parameters ($\log\xi$), distinguished by color. The `nthcomp` incident spectra are obtained from `xillvercp` by setting `refl_frac` = 0 and $kT_E = 60$keV, with $\Gamma = 1.4$ (top panel) and $\Gamma = 2.4$ (bottom panel). Other parameters are set at their default values. The incident spectrum is shown as a black dashed line for comparison.
```

## Effect of incidence angle

In `DAO`, the coronal illumination field is specified by a single incidence angle. The $I_\mathrm{inc}$ at the top boundary of the second-order radiative transfer equation is set to zero everywhere except at $\mu_{\mathrm{inc}}$.

Since the flux is defined as the first moment of the intensity:

```{math}
F_x(\tau,E) = \int_0^1 u(\tau,E,\mu) \mu d\mu
```

where $\mu = \cos\theta$. We adopt

```{math}
I = I_\mathrm{inc}\delta(\mu-\mu_\mathrm{inc})
```
The incident intensity at the top boundary can then be evaluated as:

$$
I_{\mathrm{inc}}(E) = \frac{2F_x(E)}{\mu_{\mathrm{inc}}}
$$ 

In DAO, we always calculate $F_x(E)$ first and then use above equation to derive $I_{\mathrm{inc}}$ for the boundary condition. This relation shows that the incidence angle $\mu_{\mathrm{inc}}$ directly affects the illumination intensity, which in turn influences the radiative transfer solution. Consequently, a smaller $\mu_{\mathrm{inc}}$ is expected to produce higher gas temperatures and a higher degree of ionization at the surface.


```{figure} images/dinci.svg
:width: 100%
:align: center
:name: fig-dxi

Top row: Angle-averaged emergent intensity (left) and the corresponding temperature profiles (right), color-coded by the incidence angle $\cos\theta_\mathrm{inc}$. Bottom row: Spectral ratios relative to the $\theta_{\mathrm{inc}} = 45^\circ$ case for all incidence angles. The incidence angles are sampled linearly in $\mu = \cos\theta_{\mathrm{inc}}$ between 0 and 1 using 30 bins. The ionization parameter is $\log\xi = 2.0$, and the incident radiation field (black dashed line) is a cut-off power-law with $\Gamma = 2.0$ and $E_\mathrm{cut} = 300~\mathrm{keV}$. The other parameters are set at their default values.
```
## Ions Fractions

We will run two simulations with different ionization parameters to compare the resulting ionic fractions. The primary physical parameters are set as follows:

* Ionization Parameter (log$\xi$): Three separate runs are performed, one with zeta
= 3.0 ,2.0 and 4.0.

* Incident Spectrum: The illumination is an `nthcomp` spectrum, with a photon index $\Gamma$ = 2.0 and an electron temperature $kT_e$ = 60 keV and blackbody temperture at $0.05$ keV

<details>
<summary>Click to expand plot script</summary>

```python
import matplotlib.pyplot as plt
import matplotlib as mpl 
from astropy.io import fits
import numpy as np
import pandas as pd
import matplotlib.colors as colors
from matplotlib.ticker import FixedLocator, FixedFormatter
from matplotlib.colors import Normalize,LogNorm
from matplotlib.cm import ScalarMappable
mpl.rcParams['text.usetex'] = True
plt.style.use('seaborn-v0_8-ticks')
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Computer Modern'
mpl.rcParams['font.size'] = 20
colormap = 'viridis'
def to_roman(num):
    """Converts an integer to a Roman numeral string."""
    val_syms = [
        (10, "X"), (9, "IX"), (5, "V"), (4, "IV"),
        (1, "I")
    ]
    roman_num = ""
    tens = num // 10
    roman_num += "X" * tens
    num %= 10
    for val, sym in val_syms:
        while num >= val:
            roman_num += sym
            num -= val
    return roman_num
data = pd.read_csv('./temp_a7f0a841e1_dxi5000_abund.dat',sep="\s+",header=None).iloc[-200:]
tau = np.array(pd.to_numeric(data[0],errors='coerce'))
log_tau = np.log10(tau)
midpoints = (log_tau[:-1] + log_tau[1:]) / 2.0
first_boundary = log_tau[0] - (midpoints[0] - log_tau[0])
last_boundary = log_tau[-1] + (log_tau[-1] - midpoints[-1])
boundaries_log = np.concatenate([[first_boundary], midpoints, [last_boundary]])
tau_boundaries = 10**boundaries_log

def loadabund(path):
    hdul= fits.open(path)
    abund=hdul[1].data
    return abund

def get_ionfrac(abund,ionname):
    names = [name for name in abund.columns.names if name.startswith(ionname)]
    ionfrac = np.zeros((len(names),len(tau)))
    for i,name in enumerate(names):
        ionfrac[i] = abund[name]
    ionfrac /= ionfrac.sum(axis=0)
    return ionfrac

def getinfo(abund):
    fefrac = get_ionfrac(abund,'fe_')
    ofrac = get_ionfrac(abund,'o_')
    cfrac = get_ionfrac(abund,'c_')
    ioninfo = {
        'fe':{
            "major_ions": [1, 5, 10, 15, 20, 26],
            "major_locs": [ion + 0.5 for ion in [1, 5, 10, 15, 20, 26]],
            "major_labels":[f"Fe {to_roman(i)}" for i in [1, 5, 10, 15, 20, 26]]
        },
        'o':{
            "major_ions": [1,2,3,4,5,6,7,8],
            "major_locs": [ion + 0.5 for ion in [1,2,3,4,5,6,7,8]],
            "major_labels":[f"O {to_roman(i)}" for i in [1,2,3,4,5,6,7,8]]
        },
        "c":{
            "major_ions": [1,2,3,4,5,6],
            "major_locs": [ion + 0.5 for ion in [1,2,3,4,5,6]],
            "major_labels":[f"C {to_roman(i)}" for i in [1,2,3,4,5,6]]
        }   
    }
    return fefrac,ofrac,cfrac,ioninfo


def plotcolormap(ioninfo,ax,frac,ionname):
    ion_levels = np.arange(1, frac.shape[0] + 2)

    im = ax.pcolormesh(tau_boundaries, ion_levels, frac, 
                    shading = 'flat', 
                    cmap = colormap,
                    norm = colors.LogNorm(vmin=1e-4, vmax=1.0))

    ax.set_xscale('log')
    
    major_ions = ioninfo[ionname]['major_ions']
    major_locs = ioninfo[ionname]['major_locs']
    major_labels = ioninfo[ionname]['major_labels']

    ax.yaxis.set_major_locator(FixedLocator(major_locs))
    ax.yaxis.set_major_formatter(FixedFormatter(major_labels))

    all_ions = np.arange(1, frac.shape[0] + 1)
    minor_locs = [ion + 0.5 for ion in all_ions]

    ax.yaxis.set_minor_locator(FixedLocator(minor_locs))

    ax.tick_params(axis='y', which='major', length=7, width=1.5)
    ax.tick_params(axis='y', which='minor', length=4, width=1)
    return im

fig, ax = plt.subplots(3,3,figsize=(12, 9),layout="constrained",sharex=True)
fig.supylabel(r'Ionization State')

paths = ['abund_a7f0a841e1_dxi5000_abund.fits','abund_0eae9aa8c6_dxi5000_abund.fits','abund_a2cf2f5d08_dxi5000_abund.fits']
zetas = ['4.0','3.0','2.0']

cmap = mpl.cm.get_cmap('viridis')
norm = LogNorm(vmin=1e-4, vmax=1.0)
sm = ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([]) 

for i,path in enumerate(paths):

    abund = loadabund(path)
    fefrac,ofrac,cfrac,ioninfo = getinfo(abund)

    plotcolormap(ioninfo,ax[i][0],fefrac,'fe')
    plotcolormap(ioninfo,ax[i][1],ofrac,'o')
    im=plotcolormap(ioninfo,ax[i][2],cfrac,'c')


    ax[i][0].text(0.05, 0.15, rf'log$\xi$ = {zetas[i]}', transform=ax[i][0].transAxes, verticalalignment='top', horizontalalignment='left')
    if i==len(paths)-1:
        for j in range(3):
            ax[i][j].set_xlabel(r"$\tau_{T}$")


cbar = fig.colorbar(im, ax=ax, orientation='horizontal', 
                    aspect=40, pad=0.04)
cbar.set_label("Ionization Fraction", fontsize=20)
plt.show()
```
</details>


```{figure} images/abund.svg
:width: 100%
:align: center
:name: fig-abund

Ionic fractions of iron (Fe), oxygen (O), and carbon (C) as a function of Thomson optical depth for three different ionization parameters: $\log\xi = 4.0$ (top row), $\log\xi = 3.0$ (middle row), and $\log\xi = 2.0$ (bottom row). The other parameters are the same as in Figure~\ref{fig:xivstau}.
```


## Heating and cooling rate 

we compare reflection spectra for hydrogen densities of 10$^{15}$, 10$^{16}$, and 10$^{17}$ cm$^{-3}$ under `nthcomp` illumination. Other parameters are set at their default values.

<details>
<summary>Click to expand plot script</summary>

```python
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl 
import json
import os 
from rdata import read_spec,read_spec2,read_temp,read_spec5000,readabund
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from datetime import datetime, time
import sys
def fix_scientific(x):
    if isinstance(x, str) and '-' in x and 'E' not in x and 'e' not in x:
        parts = x.split('-')
        if len(parts) == 2 and parts[0].replace('.', '').isdigit() and parts[1].isdigit():
            return float(parts[0] + 'e-' + parts[1])
    return float(x)
    
def convert_filepath(filepath):
  new_filepath = filepath.replace('spec_', 'abund_').replace('.dat', '.fits')
  return new_filepath

plt.style.use('seaborn-v0_8-ticks')
mpl.rcParams['text.usetex'] = True

mpl.rcParams['text.latex.preamble'] = r'\usepackage{mathptmx}' 
mpl.rcParams['font.family'] = 'serif'

mpl.rcParams['font.size'] = 20

ERGSEV = 1.60217662e-12 
METADATA_DIR = 'metadata'
if not os.path.isdir(METADATA_DIR):
    raise FileNotFoundError(f"元数据目录 '{METADATA_DIR}' 不存在。")
datas = os.listdir(METADATA_DIR)
def normtopflux(nh,flux):
    factor = 1e15/nh
    return flux*factor

def normhc(nh, data):
    factor = (nh/1e15)**2
    data/=factor
    return data

plot_data_list = []
for data_filename in datas:
    with open(os.path.join(METADATA_DIR, data_filename), 'r') as f:
        try:
            data = json.load(f)
            params = data['parameters']
            
            if params.get('diy_tag') == 'DNHX' and float(params['nh'])!=1e21:
                try:    
                    specfile = params['op_spec']
                    abundfile = convert_filepath(specfile)
                    tempfile = params['op_temp']

                    nh = float(params['nh'])
                    
                    # read out spec
                    ener,top,bot,inte,flux = read_spec5000(specfile)

                    # norm flux
                    flux = normtopflux(nh,flux)
                    
                    # read temperature
                    tau,temp = read_temp(tempfile)

                    # read heating-cooling
                    comh,comc,brec,toth,totc = readabund(abundfile)
                    comh = normhc(nh,comh)
                    comc = normhc(nh,comc)
                    brec = normhc(nh,brec)
                    toth = normhc(nh,toth)
                    totc = normhc(nh,totc)
   
                    task_info = {
                        'ener': ener, 
                        'eemo': flux,
                        'tau' : tau,
                        'temp': temp,
                        'top': top,
                        'nh' : float(params['nh']),
                        'compton heating':comh,
                        'compton cooling':comc,
                        'brems cooling':brec,
                        'total heating':toth,
                        'total cooling':totc
                    }
                    plot_data_list.append(task_info)
                except Exception as e:
                    print(f"fail {params['op_spec']}: {e}")
                    continue      
        except (json.JSONDecodeError, KeyError, TypeError) as e:
            print(f"Warning: skip {data_filename}: {e}")
            continue
if not plot_data_list:
    print("Error: no 'diy_tag ")
    exit()
plot_data_list.sort(key=lambda item: item['nh'])
nh_values = np.log10(np.array([item['nh'] for item in plot_data_list]))

cmap = cm.viridis   
norm = mcolors.LogNorm(vmin=nh_values.min(), vmax=nh_values.max())
model_colors = {d: mcolors.to_hex(cmap(norm(d))) for d in nh_values}

labels = {
    1e15: r'logn$_h$=15',
    1e16: r'logn$_h$=16',
    1e17: r'logn$_h$=17',
    1e18: r'logn$_h$=18',
    1e19: r'logn$_h$=19',
    1e20: r'logn$_h$=20',
    1e21: r'logn$_h$=21',
    1e22: r'logn$_h$=22'
}
fig, ax = plt.subplots(
    1, 2, 
    figsize=(3.4, 5.0),
    constrained_layout=True
)

lw = 2.0

traget_inci = np.cos(np.deg2rad(45))
for i,item in enumerate(plot_data_list):
    color = model_colors[np.log10(item['nh'])]
    label = labels[item['nh']]
    energ = item['ener']
    flux = item['eemo']
    top = item['top']
    tau=item['tau']
    temp = item['temp']

    comh = item['compton heating']
    comc = item['compton cooling']
    brec = item['brems cooling']
    toth = item['total heating']
    totc = item['total cooling']

    # if i==0:
    #     ax[0].loglog(ener, top, color='k',lw=lw,linestyle='--')

    ax[0].loglog(tau, comh , color=color, lw=1,label=label)
    ax[1].loglog(tau, comc,color=color,lw=lw,label=label)



ax[0].set_ylabel(r'Rate')
ax[0].set_xlim(0.0001,10)
ax[1].legend(loc='upper right')
ax[0].set_ylim(1e3,1e7)
ax[0].set_xticks([0.001,0.01,0.1,1,10])
ax[0].set_xticklabels(['0.001','0.01','0.1','1','10'])
ax[0].set_xlabel(r"Thomson depth $\tau_T$")
ax[0].text(0.05, 0.95, rf'Compton Heating', transform=ax[0].transAxes, verticalalignment='top', horizontalalignment='left')

ax[1].set_ylabel(r'Rate')
ax[1].set_xlim(0.0001,10)
ax[1].set_xticks([0.001,0.01,0.1,1,10])
ax[1].set_xticklabels(['0.001','0.01','0.1','1','10'])
ax[1].set_ylim(10**3,1e7)
ax[1].text(0.05, 0.95, rf'Compton Cooling', transform=ax[1].transAxes, verticalalignment='top', horizontalalignment='left')
ax[1].set_xlabel(r"Thomson depth $\tau_T$")

plt.show()
```
</details>

```{figure} images/comhc.svg
:width: 100%
:align: center
:name: fig-comhc

Compton heating and cooling rates
```

```{figure} images/PIff.svg
:width: 100%
:align: center
:name: fig-PIff

Photoionization and bremsstrahlung heating rates, along with the bremsstrahlung cooling rate.
```

## Relativistic spectrum

The `reldao` and `reldaoA` models are developed within the framework of `relxill` v2.5 {cite:p}`2014ApJ...782...76G,2025ApJ...989..168H` and `relxillA` {cite:p}`2025ApJ...989..168H`. Currently, as the full table models have not yet been generated, `reldao` serves as a prototype rather than a finalized relativistic reflection model. The code will be publicly released upon the completion of the `DAO` table models. Additionally, the relativistic spectrum is computed using the convolution model `relconv` {cite:p}`2010MNRAS.409.1534D` within `XSPEC` {cite:p}`1996ASPC..101...17A`, yielding the combined model `relconv*dao`.

The primary distinction between the convolution model `relconv*dao` and the direct implementations, `reldao` and `reldaoA`, lies in the treatment of angular dependence. The convolution model operates on the angle-averaged flux, whereas `reldao` accounts for the full angular distribution by integrating over all local emission angles at a given incidence angle. Furthermore, `relxillA` incorporates a more comprehensive treatment of angular effects. In the specific case of `reldao`, the observed flux is calculated as:

$$
F_{obs} (E_{obs}) &= \frac{1}{D^2} \int_{R_{in}}^{R_{out}} dr_e \int_0^1 g^* \frac{\pi r_e g^2}{\sqrt{g^*(1-g^*)}} \left[f^{(1)}(g^*,r_e,\theta_{obs})+f^{(2)}(g^*,r_e,\theta_{obs})\right] \\ &\times \epsilon(r_e) \langle\bar{I_e}(E_e)\rangle
$$

and for `reldaoA`, the model calculates the observed flux as:

$$
\begin{aligned}
F_{obs} (E_{obs}) &= \frac{1}{D^2} \sum_{i=0}^9 \int_{R_{in}}^{R_{out}}dr_e\int_0^1 dg^* \frac{\pi r_eg^2}{\sqrt{g^*(1-g^*)}}\left[f^{(1)}(g^*,r_e,\theta_{obs})+f^{(2)}(g^*,r_e,\theta_{obs})\right]  \\
&\times\epsilon(r_e)\bar{I}(E_e,r_e,\bar{\theta_e})\Theta(\theta_e-\theta_i)\Theta(\theta_{i+1}-\theta_e)
\end{aligned}
$$

All the physical parameters above have been explained in detail in section 2 of {cite:t}`2025ApJ...989..168H`.

```{figure} images/reldaonthcomp.svg
:name: fig-reldao
:width: 100%
:align: center

Relativistic reflection spectra calculated with `relconv*dao` (green dash-dot line), `reldao` (blue dashed line), and `reldaoA` (red solid line). The gas is illuminated by `nthcomp` with $\Gamma = 2.4$ (left column) and $\Gamma = 1.4$ (right column). Ionization parameter is set to 3.0 for both cases. The inclination angle is $\theta_{i} = 30$ deg, and all other parameters in `relxill` framework (e.g., spin, radius, etc.) are fixed at their default values.