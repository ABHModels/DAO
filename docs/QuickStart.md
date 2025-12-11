(quickstart_guide)=
# Quick Start

To run `DAO`, you must provide an input file that specifies the physical parameters and grid configurations.

We provide a template file located at `templates/input.nml.j2`.

**1. Configure the Template**

Before running a test, open the template file and change the `rfcomp_file` and `rfdataenv` parameters to match your local paths (where you downloaded the data files).

**2. Run the Simulation**

Execute the following command to run a simple test:

```bash
./dao.x templates/input.nml.j2
```

**3. Check Results**
This task will generate four output files containing the log, intensity, mean intensity, illumination, and ion abundances.

For detailed file descriptions, please refer to [Output Files](outputfile).

(runcodewithpython)=
## Run Code with Python

To streamline the workflow, we provide a Python helper script (`modelrun.py`) optimized for batch processing and parameter sweeps. This tool ensures data traceability by generating a unique hash value derived from the complete set of input parameters.

When combined with a user-defined `diy_tag` (see [Sample Output](sampleoutput)), this hash creates a unique identifier for each simulation. The script uses this identifier to:

1. Generate the model input file.

2. Archive all parameters into a JSON file for record-keeping.

This mechanism ensures that output files are strictly determined by the parameter configuration and the tag, preventing accidental overwrites. The script operates in two distinct modes:

1.  **Generation Mode**: 
    Generates the standard input file (`.nml`) and simultaneously archives the input parameters into a JSON file for record-keeping.

2.  **Slurm Submission Mode**: 
    Generates the input file and automatically submits the job to a computing cluster via the Slurm workload manager.

```{note}
**Cluster Configuration**

If you intend to use the **Slurm Submission Mode**, you must first configure the submission template. 

Please modify `templates/submission.slurm.j2` to match your specific cluster environment (e.g., partition names, time limits, and memory requirements).
```

The `MODEL_PARAMETERS` dictionary in `mopar.py` serves as the central registry for all input parameters.

```{note}
**Variable Name Mapping**

Please note that the variable names used in the Python interface differ from the internal names expected by the main Fortran program. 

To bridge this gap, the code implements an internal mapping system that automatically translates Python variables into the corresponding Fortran Namelist parameters during runtime.

```

### Bash script

The `run.sh` file serves as a practical example for executing the `modelrun.py` utility. This script is designed to be highly customizable to suit specific research needs.

* **Parameter Sweeps**: You can extend the script by adding loops to iterate over a specific range of parameter values.
* **Extended Configuration**: Users are encouraged to modify the script to include any additional parameters that are available in the `mopar.py` registry.

### Execution Script (`run.sh`)

You can use this Bash script to manage parameter sweeps or single runs.

```{code-block} bash
:caption: run.sh
:linenos:

#!/bin/bash

PYTHON_SCRIPT="modelrun.py"
DIY="User_Guide"

# --- Execution & Model Flags ---
exec=0      # 0: Generate files only, 1: Submit sbatch, 2: Run locally
debug=0     # Debug switch
model=0     # 0: Reflection, 1: Pure scattering
ity=2       # Type of angle bins and integral method

# --- Computational Grids ---
nmaxmain=15       # Max main iteration steps
nmaxrte=200       # Max radiative transfer steps
ndrt=200          # Number of depth grids (tau)
nfrt=5000         # Number of energy grids

# --- Physical Parameters ---
temp_gas=1e8      # Initial gas temperature
inci_type="nthcomp"
gamma=2           # Photon index
inmu=0.55         # Incident angle cosine
fe=0.5            # Iron abundance

# --- Loop Example ---
# You can uncomment loops to perform parameter sweeps
# for gamma in 1.0 1.2 1.4; do
    python ${PYTHON_SCRIPT} \
        --exec ${exec} \
        --debug ${debug} \
        --model ${model} \
        --ity ${ity} \
        --nmaxmain ${nmaxmain} \
        --nmaxrte ${nmaxrte} \
        --temp-gas ${temp_gas} \
        --inci-type "${inci_type}" \
        --gamma ${gamma} \
        --inmu ${inmu} \
        --fe ${fe} \
        --diy "${DIY}"
# done
```
(sampleoutput)=
### Sample Output

When running the script, you should see output indicating the configuration loading and file generation:

```{code-block} text
:caption: Terminal Output

----------------------------------------
Starting model run with DIY tag: User_Guide
----------------------------------------

--- Generated Unique Run ID: b912ee2e38_User_Guide ---

--- Final Python Configuration ---
==================================================
Current Model Configuration
==================================================
inci_file           : incident/nthcomp_g2.0_t60.txt
kernelpath          : kernel_exact5000
op_spec             : out/User_Guide/spec_b912ee2e38_User_Guide.dat
...
nmaxmain            : 15
temp_gas            : 100000000.0
fe                  : 0.5
...
h                   : 1.0
he                  : 1.0
fe                  : 0.5
ni                  : 1.0
... (other elements omitted) ...
==================================================
Parsing Fortran namelist file: templates/input.nml.j2

Successfully wrote updated parameters to 'inputFILE/inp_b912ee2e38_User_Guide.in'

--- Saving run metadata to JSON ---
Successfully saved metadata to 'metadata/meta_b912ee2e38_User_Guide.json'
----------------------------------------
Job submission script has finished.
```
