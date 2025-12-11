(Installation)=
# Installation

## Prerequisites

Before installing DAO, please ensure your system meets the following requirements:

1. **Build Tools**: `make`, `gfortran`, and `gcc`.
2. **Python**: Version 3.10 or newer.
3. **Python Libraries**: `f90nml`, `jinja2`, `hashlib`, `argparse`, and standard libraries (`os`, `subprocess`, `pathlib`).
4. **External Software**: [HEASoft](https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/) (High Energy Astrophysics Software).

## Installation Steps

1. **Clone the Repository**

   Download the source code from GitHub:

   ```bash
   git clone https://github.com/ABHModels/DAO.git
   cd DAO

2. **Prepare the Build Directory**

    Create a directory for module files:
    ```bash
    mkdir mods

3. **Configure HEASoft Paths**

    Open the `Makefile` and locate the HEADS variable. You must update this path to point to your local HEASoft installation directory.

    *Example*:

    ```bash
    HEADS = /path/to/heasoft/heasoft-6.31.1/x86_64-pclinux-gnu-libc2.17
    ```

    ```{note}
    If the libraries listed in the `LIBS` section of the `Makefile` are not found at this path, please manually verify the locations of `libcfitsio`, `libxanlib`, and `libape` on your system and update the path accordingly.

4. **Compile the Model**

    Run the make command to build the project:

    ```bash
    make

5. **Before Running: Prepare Data Files**

   The model relies on the [XSTAR](https://heasarc.gsfc.nasa.gov/docs/software/xstar/xstar.html) atomic database (specifically `adtb.fits` and `coheat.dat`) and pre-calculated Compton scattering redistribution functions.

   * **Atomic Database**: While the code can automatically locate these files if HEASoft is installed, **we strongly recommend** downloading them manually and configuring a custom path to ensure consistency. After v1.2.2, the model require a specific structure for atom database for adopt high density calculation (see [version](versionpage) for more details). 
   * **Redistribution Functions**: You must ensure the pre-calculated redistribution functions are available before running the model.

   > **Download:** You can download both the required atomic database and the Compton scattering redistribution functions (500 energy bins) from [Zenodo](https://zenodo.org/records/17845422).