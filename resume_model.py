#!/usr/bin/env python3
"""Continue a production run using its saved params, T and angular moments."""

import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def command_for(run_dir, limit, iteration=None, check=False):
    run_dir = Path(run_dir).resolve()
    if run_dir.parent != ROOT / "results":
        raise ValueError("Use a run directory directly under this repository's results/")
    with (run_dir / "params.json").open() as stream:
        params = json.load(stream)
    if params["hash"] != run_dir.name or params["test_rt"]:
        raise ValueError("The directory/hash must match a production run")
    if params.get("scattering_density") != "total_electrons_1.21_nH":
        raise ValueError("This executable uses 1.21 nH scattering targets; older models require their original executable")
    if params["ktype"] != 1 or params["maxiter"] != 500:
        raise ValueError("This executable must use the original kernel and internal RT iteration cap")
    if params.get("frac", -1) > 0 and params.get("frac_convention") != "disk_over_corona_flux_v1":
        raise ValueError("Legacy corona/disk -frac checkpoints cannot resume with the new disk/corona definition")
    # Never infer a completed iteration from thermal_status: failed runs may
    # report 101 although the latest completed history/spectrum is iteration 100.
    history = []
    for line in (run_dir / "thermal_history.dat").read_text().splitlines():
        if line.strip() and not line.lstrip().startswith("#"):
            history.append(line.split())
    if not history:
        raise ValueError("No completed outer iterations are available")
    last = int(history[-1][0])
    if iteration is None:
        iteration = last
    if iteration != last:
        raise ValueError(f"--iteration must match the last recorded iteration ({last})")
    if limit <= iteration:
        raise ValueError(f"--max-iterations must be greater than {iteration}")

    command = [str(ROOT / "maindaocl")]
    if params.get("thermal_update") != "full" or not params.get("log_depth", False):
        raise ValueError("This executable resumes only full-temperature, logarithmic-grid checkpoints")
    for key, flag in (
        ("corona", "-corona"), ("nh", "-nh"), ("zeta", "-zeta"),
        ("frac", "-frac"), ("incidence", "-incidence"), ("Afe", "-Afe"),
        ("Gamma", "-Gamma"), ("E_cut", "-Ecut"), ("E_lo_cut", "-E_low_cut"),
        ("kT_e", "-kT_e"), ("kT_bb", "-kT_bb"), ("taup", "-taup"),
        ("kT_disk", "-kT_disk"), ("tau_slab", "-tau"),
    ):
        command.extend([flag, str(params[key])])
    for key, flag in (("angsca", "-angsca"),
                      ("save_oxygen", "-O"), ("save_iron", "-Fe")):
        command.extend([flag, str(int(params.get(key, False)))])
    if params.get("label"):
        command.extend(["-label", params["label"]])
    command.extend(["-resume_iter", str(iteration), "-max_outer_iter", str(limit)])
    if check:
        command.append("-resume_check")
    return command, params, iteration


def runtime_environment():
    env = os.environ.copy()
    cache = ROOT.parent / "DAO_KERNEL_TABLE"
    if cache.is_dir():
        env.setdefault("COMPTON_CACHE_DIR", str(cache))
    # Match this checkout's existing HEASoft build; inherit any user overrides.
    xspec_lib = Path("/Users/xe26734/Downloads/heasoft-6.37.1/aarch64-apple-darwin25.6.0/lib")
    if sys.platform == "darwin" and xspec_lib.is_dir():
        env["DYLD_LIBRARY_PATH"] = str(xspec_lib) + ":/opt/homebrew/lib" + (
            ":" + env["DYLD_LIBRARY_PATH"] if env.get("DYLD_LIBRARY_PATH") else "")
    return env


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", help="Existing results/<hash> directory")
    parser.add_argument("--max-iterations", type=int, default=300,
                        help="Total outer-iteration cap, including saved iterations (default: 300)")
    parser.add_argument("--iteration", type=int, help="Saved iteration; defaults to the latest history row")
    parser.add_argument("--check", action="store_true", help="Validate saved state without calculating or writing output")
    parser.add_argument("--dry-run", action="store_true", help="Print the command without running it")
    args = parser.parse_args()
    try:
        command, params, iteration = command_for(args.run_dir, args.max_iterations, args.iteration, args.check)
        print(f"{params['hash']}: resume {iteration} → cap {args.max_iterations} "
              f"({'angdep' if params['angsca'] else 'avgkernel'}; mu_inc={params['incidence']})", flush=True)
        if args.dry_run:
            print(shlex.join(command))
            return 0
        if not Path(command[0]).is_file():
            raise ValueError("Build the updated executable first: make -j4")
        # The binary also validates all physics metadata and full T/J grids.
        if args.check:
            return subprocess.run(command, cwd=ROOT, env=runtime_environment()).returncode
        preflight = subprocess.run(command + ["-resume_check"], cwd=ROOT, env=runtime_environment())
        if preflight.returncode:
            return preflight.returncode
        os.chdir(ROOT)
        os.execve(command[0], command, runtime_environment())
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        print(f"Cannot resume: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
