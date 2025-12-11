# mopar/__init__.py

# ==============================================================================
#                 MODEL PARAMETERS (Single Source of Truth)
# ==============================================================================
# Author: Yimin Huang
# Affiliation: Fudan University
#
#
# This list defines every parameter used by the model.
# It is used to automatically generate command-line arguments,
# sync with the Fortran namelist, and generate unique run IDs.
# ==============================================================================
import os
MODEL_PARAMETERS = [
    # --- I/O & File Path Parameters ---
    {
        'py_name': 'inci_file',     'type': str,   'default': 'incident/nthcomp_g2.0_t60.txt',
        'help': 'Incident spectrum file name [Necessary for type of incident is "file"]'
    },
    {
        'py_name': 'kernelpath',    'type': str,   'default': 'kernel_exact5000',
        'help': 'Compton scattering Redistribution function file, the energy bins need to be as same in model'
    },
    {
        'py_name': 'dataenv',       'type': str,   'default': 'data',
        'help': 'Atomic database and Compton heating-cooling file'
    },
    {
        'py_name': 'op_spec',       'type': str,   'default': "spec.dat",
        'help': 'output spectrum file'
    },
    {
        'py_name': 'op_temp',       'type': str,   'default': "temp.dat",
        'help': 'output temperature file'
    },
    {
        'py_name': 'op_inte',       'type': str,   'default': "inte.dat",
        'help': 'output intensity file'
    },
    {
        'py_name': 'op_log',        'type': str,   'default': "model.log",
        'help': 'output log file'
    },

    {
        'py_name': 'op_emis',        'type': str,   'default': "emis.dat",
        'help': 'output emissivity file'
    },

    {
        'py_name': 'op_abund',        'type': str,   'default': "abud.fits",
        'help': 'output abund file'
    },

    # --- Control & Numerical Precision Parameters ---
    {
        'py_name': 'nmaxmain',      'type': int,   'default': 15,
        'help': 'Number of maximum main iteration steps'
    },
    {
        'py_name': 'nmaxrte',       'type': int,   'default': 100,
        'help': 'Number of maximum radiative transfer iteration steps'
    },
    {
        'py_name': 'ndrt',          'type': int,   'default': 200,
        'help': 'Number of depth grids'
    },
    {
        'py_name': 'nfrt',          'type': int,   'default': 5000,
        'help': 'Number of energy grids'
    },
    {
        'py_name': 'nmurt',         'type': int,   'default': 10,
        'help': 'Number of angle grids'
    },
    {
        'py_name': 'ppemin',        'type': float, 'default': 0.1,
        'help': 'Minimum energy in eV'
    },
    {
        'py_name': 'ppemax',        'type': float, 'default': 900000.0,
        'help': 'Before maximum energy in eV'
    },
    {
        'py_name': 'ppemax2',       'type': float, 'default': 1000000.0,
        'help': 'Maximum energy in eV'
    },
    {
        'py_name': 'taumin',        'type': float, 'default': 0.0001,
        'help': 'Minimum Thomson optical depth'
    },
    {
        'py_name': 'taumaxrt',      'type': float, 'default': 10.0,
        'help': 'Maximum Thomson optical depth'
    },
    {
        'py_name': 'mumin',         'type': float, 'default': 0.0,
        'help': 'Minimum cosine value'
    },
    {
        'py_name': 'mumax',         'type': float, 'default': 1.0,
        'help': 'Maximum cosine value'
    },

    # --- Model Flags & Physical Parameters ---
    {
        'py_name': 'ity',           'type': int,   'default': 2,
        'help': 'Type of angle bins and integral method for angle'
    },
    {
        'py_name': 'debug',         'type': int,   'default': 0,
        'help': 'Debug switch'
    },
    {
        'py_name': 'it_temp',       'type': int,   'default': 99,
        'help': 'Maximum iteration steps for thermal equilibrium'
    },
    {
        'py_name': 'temp_gas',      'type': float, 'default': 1e8,
        'help': 'Initial gas temperature'
    },
    {
        'py_name': 'temp_gas_unit', 'type': str,   'default': 'k',
        'help': 'Unit of initial gas temperature'
    },
    {
        'py_name': 'nh',            'type': float, 'default': 1.0e15,
        'help': 'Hydrogen density'
    },
    {
        'py_name': 'model',         'type': int,   'default': 0,
        'help': 'Model selector [1: Pure scattering 0: reflecion]'
    },

    # --- Incident Spectrum Parameters ---
    {
        'py_name': 'inci_type',     'type': str,   'default': 'powerlaw',
        'help': 'Type of corona illumination'
    },
    {
        'py_name': 'gamma',         'type': float, 'default': 2.0,
        'help': 'Photon index or Blackbody temperature in eV'
    },
    {
        'py_name': 'hcut',          'type': float, 'default': 300E3,
        'help': 'High energy cut-off for cut-off power law illumination or electrons temperture'
    },
    {
        'py_name': 'ktbb_nth',      'type': float,  'default':0.1e3,
        'help': 'Blackbody temperature when incident type is Comptonization source'
    },
    {
        'py_name': 'inmu',          'type': float, 'default': 0.7,
        'help': 'Incidence angle'
    },
    {
        'py_name': 'sbot',          'type': int,   'default': 0,
        'help': 'Switch of bottom illumination [Thermal disk radiation]'
    },
    {
        'py_name': 'ktbb',          'type': float, 'default': 350.0,
        'help': 'Temperature of disk radiation'
    },
    {
        'py_name': 'fx_frac',       'type': float, 'default': 1.0,
        'help': 'Flux of top illumiantion over total illumination [defined by ionizaiton parameter and hydrogen density]'
    },
    {
        'py_name': 'zeta',          'type': float, 'default': 3.0,
        'help': 'Log10 of ionization parameter'
    },

    # --- Elemental Abundances ---
    # Defines the abundance of each element relative to Hydrogen.
    { 'py_name': 'h',   'type': float, 'default': 1.0, 'help': 'Abundance of Hydrogen. This is the reference element.' },
    { 'py_name': 'he',  'type': float, 'default': 1.0, 'help': 'Abundance of Helium relative to Hydrogen.' },
    { 'py_name': 'li',  'type': float, 'default': 0.0, 'help': 'Abundance of Lithium relative to Hydrogen.' },
    { 'py_name': 'be',  'type': float, 'default': 0.0, 'help': 'Abundance of Beryllium relative to Hydrogen.' },
    { 'py_name': 'ba',  'type': float, 'default': 0.0, 'help': 'Abundance of Barium relative to Hydrogen.' },
    { 'py_name': 'c',   'type': float, 'default': 1.0, 'help': 'Abundance of Carbon relative to Hydrogen.' },
    { 'py_name': 'n',   'type': float, 'default': 1.0, 'help': 'Abundance of Nitrogen relative to Hydrogen.' },
    { 'py_name': 'o',   'type': float, 'default': 1.0, 'help': 'Abundance of Oxygen relative to Hydrogen.' },
    { 'py_name': 'f',   'type': float, 'default': 0.0, 'help': 'Abundance of Fluorine relative to Hydrogen.' },
    { 'py_name': 'ne',  'type': float, 'default': 1.0, 'help': 'Abundance of Neon relative to Hydrogen.' },
    { 'py_name': 'na',  'type': float, 'default': 0.0, 'help': 'Abundance of Sodium relative to Hydrogen.' },
    { 'py_name': 'mg',  'type': float, 'default': 1.0, 'help': 'Abundance of Magnesium relative to Hydrogen.' },
    { 'py_name': 'al',  'type': float, 'default': 0.0, 'help': 'Abundance of Aluminum relative to Hydrogen.' },
    { 'py_name': 'si',  'type': float, 'default': 1.0, 'help': 'Abundance of Silicon relative to Hydrogen.' },
    { 'py_name': 'p',   'type': float, 'default': 0.0, 'help': 'Abundance of Phosphorus relative to Hydrogen.' },
    { 'py_name': 's',   'type': float, 'default': 1.0, 'help': 'Abundance of Sulfur relative to Hydrogen.' },
    { 'py_name': 'cl',  'type': float, 'default': 0.0, 'help': 'Abundance of Chlorine relative to Hydrogen.' },
    { 'py_name': 'ar',  'type': float, 'default': 1.0, 'help': 'Abundance of Argon relative to Hydrogen.' },
    { 'py_name': 'k',   'type': float, 'default': 0.0, 'help': 'Abundance of Potassium relative to Hydrogen.' },
    { 'py_name': 'ca',  'type': float, 'default': 1.0, 'help': 'Abundance of Calcium relative to Hydrogen.' },
    { 'py_name': 'sca', 'type': float, 'default': 0.0, 'help': 'Abundance of Scandium relative to Hydrogen.' },
    { 'py_name': 'ti',  'type': float, 'default': 0.0, 'help': 'Abundance of Titanium relative to Hydrogen.' },
    { 'py_name': 'va',  'type': float, 'default': 0.0, 'help': 'Abundance of Vanadium relative to Hydrogen.' },
    { 'py_name': 'cr',  'type': float, 'default': 0.0, 'help': 'Abundance of Chromium relative to Hydrogen.' },
    { 'py_name': 'mn',  'type': float, 'default': 0.0, 'help': 'Abundance of Manganese relative to Hydrogen.' },
    { 'py_name': 'fe',  'type': float, 'default': 1.0, 'help': 'Abundance of Iron relative to Hydrogen.' },
    { 'py_name': 'co',  'type': float, 'default': 0.0, 'help': 'Abundance of Cobalt relative to Hydrogen.' },
    { 'py_name': 'ni',  'type': float, 'default': 1.0, 'help': 'Abundance of Nickel relative to Hydrogen.' },
    { 'py_name': 'cu',  'type': float, 'default': 0.0, 'help': 'Abundance of Copper relative to Hydrogen.' },
    { 'py_name': 'zn',  'type': float, 'default': 0.0, 'help': 'Abundance of Zinc relative to Hydrogen.' },
]

# ==============================================================================
#                 2. Python Name -> Fortran Namelist Name MAP
# ==============================================================================
# This dictionary acts as a translator between the Python script's
# parameter names and the actual variable names inside the Fortran
# namelist file.
# ==============================================================================

PARAM_MAP = {
    # --- I/O & File Path Parameters ---
    'inci_file':     'rfinci_file',
    'kernelpath':    'rfcomp_file',  # Special mapping
    'dataenv':       'rfdataenv',
    'op_spec':       'rfop_spec',
    'op_temp':       'rfop_temp',
    'op_inte':       'rfop_inte',
    'op_log':        'rfop_log',
    'op_emis':       'rfop_emis',
    'op_abund':      'rfabdfits',         

    # --- Control & Numerical Precision Parameters ---
    'nmaxmain':      'rfnmaxmain',
    'nmaxrte':       'rfnmaxrte',
    'ndrt':          'ndrt',
    'nfrt':          'nfrt',
    'nmurt':         'nmurt',
    'ppemin':        'rfppemin',
    'ppemax':        'rfppemax',
    'ppemax2':       'rfppemax2',
    'taumin':        'rftaumin',
    'taumaxrt':      'rftaumaxrt',
    'mumin':         'rfmumin',
    'mumax':         'rfmumax',

    # --- Model Flags & Physical Parameters ---
    'ity':           'ity',
    'debug':         'rfdebug',
    'it_temp':       'rfit_temp',
    'temp_gas':      'rftemp_gas',
    'temp_gas_unit': 'rftemp_gas_unit',
    'nh':            'rfnh',
    'model':         'rfmodel',

    # --- Incident Spectrum Parameters ---
    'inci_type':     'rfstinci',  # Special mapping
    'gamma':         'rfgamma',
    'hcut':          'rfhcut',
    'inmu':          'rfinmu',
    'sbot':          'rfsbot',
    'ktbb':          'rfktbb',
    'fx_frac':       'rffx_frac',
    'zeta':          'rfzeta',
    "ktbb_nth":      'ktb_nthcomp',

    # --- Elemental Abundances ---
    'h':   'rfh',    'he':  'rfhe',   'li':  'rfli',   'be':  'rfbe',
    'ba':  'rfba',   'c':   'rfc',    'n':   'rfn',    'o':   'rfo',
    'f':   'rff',    'ne':  'rfne',   'na':  'rfna',   'mg':  'rfmg',
    'al':  'rfal',   'si':  'rfsi',   'p':   'rfp',    's':   'rfs',
    'cl':  'rfcl',   'ar':  'rfar',   'k':   'rfk',    'ca':  'rfca',
    'sca': 'rfsca',  'ti':  'rfti',   'va':  'rfva',   'cr':  'rfcr',
    'mn':  'rfmn',   'fe':  'rffe',   'co':  'rfco',   'ni':  'rfni',
    'cu':  'rfcu',   'zn':  'rfzn',
}
TEMPLATE_DIR = 'templates'

INPUT_NML_FILE = os.path.join(TEMPLATE_DIR, 'input.nml.j2')
INPUT_SLURM_FILE = os.path.join(TEMPLATE_DIR, 'submission.slurm.j2')
