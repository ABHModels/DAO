# ==============================================================================
#                               Model Runner
# ==============================================================================
#
# Author: Yimin Huang
# Affiliation: Fudan University
#
# Description:
#   This Python script serves as a powerful wrapper for a Fortran-based scientific
#   model. It automates the entire workflow of running the model, from parameter
#   configuration to job submission and metadata tracking. The script is designed
#   to enhance reproducibility, scalability, and ease of use for conducting
#   numerical experiments.
#
# Key Features:
#   1.  Dynamic Configuration: Parses command-line arguments to configure model
#       parameters, eliminating the need for manual editing of input files.
#   2.  Automatic File Generation: Creates unique, hash-based filenames for each
#       run to prevent conflicts and ensure traceability. It generates the necessary
#       Fortran namelist input files and, optionally, Slurm submission scripts
#       from templates.
#   3.  Flexible Execution Modes: Supports multiple execution modes:
#       - Generate input files only (for manual runs or inspection).
#       - Submit the job to a Slurm cluster via `sbatch`.
#       - Run the model directly in the local terminal.
#   4.  Metadata Management: For each run, it records all input parameters,
#       output file paths, and run information into a structured JSON file,
#       providing a clear and queryable record of every experiment.
#   5.  Dependency Validation: Performs basic checks to ensure that required
#       parameters are provided based on the chosen model options (e.g.,
#       requiring a filename when the input type is 'file').
#
# Workflow:
#   1.  Parse command-line arguments to get model parameters.
#   2.  Generate a unique run ID (hash) based on the input parameters.
#   3.  Create all necessary output directories.
#   4.  Generate and write the Fortran namelist input file (`.in`).
#   5.  Based on the `--exec` flag, either:
#       a. Do nothing further.
#       b. Create a Slurm script (`.sh`) and submit it.
#       c. Execute the model binary (`.x`) directly.
#   6.  Save a comprehensive metadata JSON file for the run.
#
# ==============================================================================
import f90nml
import os
import argparse
import subprocess
from pathlib import Path
from jinja2 import Environment, FileSystemLoader
from mopar import MODEL_PARAMETERS,PARAM_MAP,INPUT_NML_FILE,INPUT_SLURM_FILE
import json  
from datetime import datetime 
import hashlib 

# ==============================================================================
#                 Utility Classes and Functions
# ==============================================================================
class ModelConfig:
    """
    A class to manage and store all model configuration parameters.
    It is initialized directly from parsed command-line arguments.
    """
    def __init__(self, args):
        """
        Initializes the configuration from parsed argparse arguments.
        
        Args:
            args (argparse.Namespace): The result of parser.parse_args().
        """
        for key, value in vars(args).items():
            setattr(self, key, value)

    def update(self, **kwargs):
        """Updates the configuration using keyword arguments."""
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                print(f"Warning: '{key}' is not a valid configuration parameter.")

    def __str__(self):
        """
        Provides a formatted string representation of the current configuration,
        maintaining the logical order from MODEL_PARAMETERS.
        """
        lines = ["="*50, "Current Model Configuration", "="*50]
        
        # 直接遍历 MODEL_PARAMETERS 列表，这样就能保持我们想要的顺序
        for param_info in MODEL_PARAMETERS:
            py_name = param_info['py_name']
            
            # 从 self (也就是 ModelConfig 对象) 中获取对应参数的值
            if hasattr(self, py_name):
                value = getattr(self, py_name)
                
                # 只打印那些有值的参数，避免打印出 None
                if value is not None:
                    lines.append(f"{py_name:<20}: {value}")
        
        lines.append("="*50)
        return "\n".join(lines)

def read_namelist(filename):
    """Reads a Fortran namelist file and returns a Namelist object."""
    print(f"Parsing Fortran namelist file: {filename}")
    try:
        return f90nml.read(filename)
    except FileNotFoundError:
        print(f"Error: Namelist file '{filename}' not found.")
        return None

def sync_config_to_nml(config, nml, nml_group='input_params'):
    """Syncs parameters from a ModelConfig object to an f90nml Namelist object."""
    if nml_group not in nml:
        print(f"Warning: Namelist group '{nml_group}' not found in the file.")
        return

    params_to_update = nml[nml_group]
    for py_name, py_value in config.__dict__.items():
        f90_name = PARAM_MAP.get(py_name)
        if f90_name and f90_name in params_to_update:
            if py_value is not None and params_to_update[f90_name] != py_value:
                # print(f"Updating '{f90_name}': {params_to_update[f90_name]} -> {py_value} (from Python '{py_name}')")
                params_to_update[f90_name] = py_value

def create_slurm_script(template_path, output_path, job_name, nml_input_file):
    """Creates a specific Slurm submission script from a template."""
    print(f"\n--- Creating Slurm script for job: {job_name} ---")
    template_dir = os.path.dirname(template_path) or '.'
    env = Environment(loader=FileSystemLoader(template_dir), trim_blocks=True, lstrip_blocks=True)
    slurm_context = {'job_name': job_name, 'nml_input_file': nml_input_file}

    try:
        template = env.get_template(os.path.basename(template_path))
        rendered_script = template.render(slurm_context)
        with open(output_path, 'w') as f:
            f.write(rendered_script)
        print(f"Successfully created Slurm submission script: '{output_path}'")
        print(f"  -> Job Name: {job_name}")
        print(f"  -> Input File: {nml_input_file}")
    except Exception as e:
        print(f"Error creating Slurm script: {e}")
def save_metadata_to_json(filepath, metadata):
    """
    Save run metadata to a JSON file.

    Args:
        filepath (str): Path to the output JSON file.
        metadata (dict): Metadata to save.
    """
    print(f"\n--- Saving run metadata to JSON ---")
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    try:
        with open(filepath, 'w') as f:
            json.dump(metadata, f, indent=4)
        print(f"Successfully saved metadata to '{filepath}'")
    except Exception as e:
        print(f"Error saving metadata to JSON: {e}")
# ==============================================================================
#                 Main Workflow
# ==============================================================================
def main(args):
    """
    Main workflow to generate Fortran input files and execute the model.
    """
    run_start_time = datetime.now().isoformat()

    # 1. Initialize config object from command-line arguments
    config = ModelConfig(args)

    # 2. Validate parameter dependencies
    if config.inci_type == 'file' and config.inci_file is None:
        parser.error('When --inci-type is "file", the --inci-file argument is required.')
    if config.inci_type == 'powlaw' and config.gamma is None:
        parser.error('When --inci-type is "powlaw", the --gamma argument is required.')
    if config.inci_type == 'cutoff' and (config.gamma is None or config.hcut is None):
        parser.error('When --inci-type is "cutoff", the --gamma and --hcut arguments are required.')
    if config.kernelpath is None:
        parser.error('The --kernelpath argument is required. Please specify the path to the kernel data file.')


    # 3. Dynamically generate a base filename from key parameters
    params_for_hashing = {}
    for param_info in MODEL_PARAMETERS:
        py_name = param_info['py_name']
        if not py_name.startswith('op_'):
            params_for_hashing[py_name] = getattr(config, py_name)
    canonical_string = json.dumps(params_for_hashing, sort_keys=True)
    hasher = hashlib.sha1(canonical_string.encode('utf-8'))
    unique_hash = hasher.hexdigest()[:10] # e.g., 'a9b1c8d0e7'
    endwith = f'{unique_hash}'
    if args.diy:
        endwith += f'_{args.diy}'
    print(f"\n--- Generated Unique Run ID: {endwith} ---")

    if not os.path.exists('out'):
        os.makedirs('out')
    
    if not os.path.exists('debug'):
        os.makedirs('debug')

    dirname = os.path.join("out",args.diy)
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    
    # 4. Update the config object with the generated output filenames
    config.update(
        op_spec=os.path.join(dirname,f'spec_{endwith}.dat'),
        op_temp=os.path.join(dirname,f'temp_{endwith}.dat'),
        op_inte=os.path.join(dirname,f'inte_{endwith}.dat'),
        op_log=os.path.join(dirname,f'model_{endwith}.log'),
        op_emis=os.path.join(dirname,f'emis_{endwith}.dat'),
        op_abund=os.path.join(dirname,f'abund_{endwith}.fits')
    )
    print("\n--- Final Python Configuration ---")
    print(config)

    # 5. Create and write the Fortran Namelist file
    if not os.path.exists('inputFILE'):
        os.makedirs('inputFILE')
    OUTPUT_NML_FILE = os.path.join('inputFILE', f'inp_{endwith}.in')
    nml = read_namelist(INPUT_NML_FILE)
    if nml:
        sync_config_to_nml(config, nml, nml_group='input_params')
        nml.write(OUTPUT_NML_FILE, force=True)
        print(f"\nSuccessfully wrote updated parameters to '{OUTPUT_NML_FILE}'")

    # 6. Execute or submit the job based on the --exec argument
    if args.exec == 1: # Submit with sbatch
        if not os.path.exists('slurm_script'):
            os.makedirs('slurm_script')
        OUTPUT_SLURM_FILE = os.path.join('slurm_script', f'sl_{endwith}.sh')
        job_name = f'{endwith}'
        create_slurm_script(INPUT_SLURM_FILE, OUTPUT_SLURM_FILE, job_name, OUTPUT_NML_FILE)
        print(f"\nSubmitting job with command: sbatch {OUTPUT_SLURM_FILE}")
        try:
            result = subprocess.run(['sbatch', OUTPUT_SLURM_FILE], check=True, capture_output=True, text=True)
            print("Job submitted successfully.")
            print("STDOUT:", result.stdout)
        except subprocess.CalledProcessError as e:
            print(f"Error submitting job with sbatch.", e)
            print("STDERR:", e.stderr)
    elif args.exec == 2: # Run locally
        print(f"\nRunning model in terminal...")
        try:
            # Assumes the executable is in the current directory
            result = subprocess.run(['./dao.x', OUTPUT_NML_FILE], check=True, capture_output=True, text=True)
            print("--- [SUCCESS] ---")
            print("STDOUT:", result.stdout)
        except subprocess.CalledProcessError as e:
            print("--- [ERROR] ---")
            print("\n--- STDOUT from failed process ---\n", e.stdout)
            print("\n--- STDERR from failed process ---\n", e.stderr)

    config_dict = {p['py_name']: getattr(config, p['py_name']) for p in MODEL_PARAMETERS}
    config_dict['exec_mode'] = args.exec
    config_dict['diy_tag'] = args.diy
    run_metadata = {
        'run_info': {
            'start_time_iso': run_start_time,
            'status': 'submitted' if args.exec in [1, 2] else 'generated_only',
            'run_id': endwith 
        },
        'parameters': config_dict,
        'output_files': {
            'namelist_input': OUTPUT_NML_FILE,
            'spectrum': config.op_spec,
            'temperature': config.op_temp,
            'intensity': config.op_inte,
            'log': config.op_log,
        }
    }

    if not os.path.exists('metadata'):
        os.makedirs('metadata')
    OUTPUT_JSON_FILE = os.path.join('metadata', f'meta_{endwith}.json')
    
    run_metadata['output_files']['metadata_json'] = OUTPUT_JSON_FILE

    save_metadata_to_json(OUTPUT_JSON_FILE, run_metadata)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Generate a Fortran namelist file and run the model.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter # Automatically show default values
    )
    # --- Automatically create command-line arguments from MODEL_PARAMETERS ---
    for param in MODEL_PARAMETERS:
        # Convert Python-style 'py_name' to CLI-style '--py-name'
        cli_name = '--' + param['py_name'].replace('_', '-')
        parser.add_argument(cli_name, type=param['type'], default=param['default'], help=param['help'])

    # --- Add script-specific arguments (not part of the model config) ---
    parser.add_argument('--exec', type=int, default=0, choices=[0, 1, 2],
                        help='Execution mode. 0: Generate files only. 1: Submit with sbatch. 2: Run locally.')
    parser.add_argument('--diy', type=str, default=None,
                        help='A custom string to append to output filenames for easy identification.')
    
    args = parser.parse_args()
    main(args)