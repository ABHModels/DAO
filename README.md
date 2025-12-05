# DAO: An Open-Source X-ray Reflection Model

**DAO** is a comprehensive, open-source model for simulating X-ray reflection spectra from photoionized accretion disks.
It solves the radiative transfer equation using the robust Feautrier method and incorporates atomic data from the powerful XSTAR database.

For a detailed description of the underlying physics and methodology, please refer to our accompanying paper.

---

##  Clone the Repository

You can clone the latest version of DAO using:

```bash
git clone https://github.com/ABHModels/DAO.git
cd DAO
```


## Prerequisites

To compile and run DAO, you will need the following software installed on your system:

*   **A Fortran Compiler**: A modern Fortran compiler, such as `gfortran`.
*   **Python 3**: Along with standard libraries.
*   **HEAsoft**: A functioning installation of the [HEAsoft software suite](https://heasarc.gsfc.nasa.gov/lheasoft/). The model relies on its libraries for FITS file handling and other functionalities.

## Installation

The model uses a `Makefile` to simplify the compilation process. Please follow these steps:

1.  **Configure the Makefile:**
    Open the `Makefile` in a text editor. You must locate and update the variable that specifies the path to your HEAsoft installation. Ensure all required libraries (`LIBS`) listed in the Makefile are available in your environment.

2.  **Compile the Code:**
    From the root directory of the project, run the `make` command:
    ```bash
    mkdir mods
    make
    ```
    If the compilation is successful, an executable file (e.g., `dao.x`) will be created in the main directory.

## Running the Model

We provide a powerful Python script (`modelrun.py`) to help you configure parameters, manage runs, and explore the parameter space efficiently.

### Before You Run: Essential Data Files

Before launching a simulation, please ensure the following data files are correctly prepared and accessible:

1.  **Atomic Database (`atdb.fits`)**: The model requires an atomic database in FITS format. ***HEAsoft*** already includes it, but make sure the atdb.fits file is located in the directory specified by your input parameters if you want to use differenet version from ***HEAsoft***.

2.  **Compton Scattering Redistribution File**: DAO requires a pre-calculated file containing the Compton scattering redistribution matrix. You must generate this file beforehand and provide the correct path to it in the model's parameters.

**Consistent Energy Bins**: **This is critical.** The energy grid used to generate the Compton scattering file **must** be identical to the energy grid used in the model simulation. We provide a utility in the `get_redistribution/` directory to help you calculate redistribution function with a consistent grid structure. Please ensure the number of energy bins is the same for both the redistribution file and your model run.

### Executing a Run

You can execute the model using the provided Python wrapper script. See [User guide](DAO_User_Guide.pdf), Chapter 3.3 for more details. 

## Physics behind DAO

We write a [User guide](DAO_User_Guide.pdf), which try to include all of the model. Please see it for more details and examples.


## LICENSE
This project is released under the MIT License.
See the [LICENSE](LICENSE) file for details.
