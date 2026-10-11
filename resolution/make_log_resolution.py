#!/usr/bin/env python3
"""Write a Cloudy continuum mesh configuration with constant log spacing."""

import argparse
import math
from pathlib import Path

RYD_EV = 13.605693


def positive_float(value):
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("must be finite and positive")
    return number


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolution", type=positive_float, default=300.0,
                        help="Cloudy resolving-power parameter R (default: 300)")
    parser.add_argument("--emin-ev", type=positive_float, default=1.0,
                        help="lower energy bound in eV (default: 1)")
    parser.add_argument("--emax-kev", type=positive_float, default=1000.0,
                        help="upper energy bound in keV (default: 1000)")
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).resolve().with_name("log_resolution.ini"),
                        help="output path (default: next to this script)")
    args = parser.parse_args()
    upper_ev = args.emax_kev * 1000.0
    if not math.isfinite(upper_ev) or upper_ev <= args.emin_ev:
        parser.error("upper energy must be finite and greater than lower energy")

    # Each row gives an UPPER bound and the resolution BELOW that bound.
    # Keep the existing smooth_hump.ini settings outside the selected range.
    args.output.write_text(
        "10 08 08 //the magic number for this format file\n"
        "# Constant logarithmic spacing within the selected energy range.\n"
        f"# E_min={args.emin_ev:g} eV  E_max={args.emax_kev:g} keV\n"
        f"# R={args.resolution:g}; nominal delta(ln E)=1/R.\n"
        "# Format: upper_limit_Ryd  resolving_power\n"
        f"{args.emin_ev / RYD_EV:.12g}  10\n"
        f"{upper_ev / RYD_EV:.12g}  {args.resolution:.12g}\n"
        "0  33.333333\n",
        encoding="utf-8",
    )
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
