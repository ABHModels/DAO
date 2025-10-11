#!/bin/bash

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                             #
#  Description: This is a runner script for modelrun.py, containing all       #
#               the input parameters for the model.                           #
#                                                                             #
#  Author: Yimin Huang, Fudan University (Updated based on mopar config)      #
#                                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#=======================================================================
# Script to run the model with a comprehensive set of parameters.
# Variable names are synchronized with the 'mopar' Python module.
#=======================================================================
PYTHON_SCRIPT="modelrun.py" # <-- Make sure this is your main Python script name
DIY="tbnth"

# --- Execution & Model Flags ---
exec=1      # (0: Generate files only, 1: Submit with sbatch, 2: Run locally)
debug=0     # Debug switch
model=0     # (0: Reflection, 1: Pure scattering)
ity=2       # Type of angle bins and integral method

# --- Computational Grids ---
nmaxmain=15        # Number of maximum main iteration steps
nmaxrte=200        # Number of maximum radiative transfer iteration steps
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
inci_type="nthcomp"               # Type of corona illumination ('powerlaw', 'file', 'cutoff')
gamma=2                           # Photon index or Blackbody temperature in eV
hcut=60                           # High energy cut-off for cut-off power law illumination [eV] or plasam temperature
inmu=0.55                         # Cosine of the incident illumination angle
inci_file="incident/nthcomp_g2.0_t60.txt" # Path to incident spectrum file

# --- Bottom Illumination (Disk) ---
sbot=0              # Switch for bottom illumination (0: off, 1: on)
ktbb=350.0          # Temperature of the bottom blackbody source [eV]
fx_frac=1e-10       # Flux ratio of top illumination to total illumination

# --- Atmosphere Properties ---
zeta=3.0            # Log10 of ionization parameter (xi)
nh=1e15             # Hydrogen number density [cm^-3]
fe=0.5              # Iron abundance relative to Hydrogen

# --- Paths ---
kernelpath="kernel_exact5000" # Path to the Compton scattering redistribution file

#=======================================================================
# --- Execute the Model ---
# The script will now call the Python program with all parameters defined above.
# The parameter names (--parameter-name) match the 'py_name' in mopar.
#=======================================================================

echo "----------------------------------------"
echo "Starting model run with DIY tag: ${DIY}"
echo "----------------------------------------"


# Example loop: You can loop over any parameter.
# Here, we loop over the incident angle 'inmu'.
inmu=0.55
zeta=3.0
fe=0.5 
for gamma in 1.0 1.2 1.4 1.6 1.8
do
for hcut in 1.0 10 20 50 80 100 150 200 250 300 350 400 
do    
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
echo "Job submission script has finished."
echo "----------------------------------------"