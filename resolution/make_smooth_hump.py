# Smooth-hump resolving-power generator for the Compton RT code.
#
# Author: Yimin Huang (Fudan University; University of Bristol)
# Email:  huangym23@m.fudan.edu.cn

import numpy as np

# Parameters
R_base   = 100.0
R_peak   = 1000.0
E_center = 50.0      # keV
sigma    = 0.5       # dex
N_zones  = 80        # number of log-spaced intervals -> 81 grid points
E_lo_eV  = 1.0       # 1 eV
E_hi_keV = 1000.0    # 1000 keV
Ryd_eV   = 13.605693  # 1 Ryd in eV
closing_R = 33.333333
first_R   = 10.0      # R at the lowest edge

# Energy grid in Ryd (log-spaced)
E_lo_Ryd = E_lo_eV / Ryd_eV
E_hi_Ryd = (E_hi_keV * 1000.0) / Ryd_eV
E_Ryd = np.logspace(np.log10(E_lo_Ryd), np.log10(E_hi_Ryd), N_zones + 1)

# Convert to keV for the resolution formula
E_keV = E_Ryd * Ryd_eV / 1000.0

# Smooth hump resolving power
z = (np.log10(E_keV) - np.log10(E_center)) / sigma
R = R_base + (R_peak - R_base) * np.exp(-0.5 * z * z)

# Override first entry with the low-edge value used in the file
R[0] = first_R

with open("smooth_hump.ini", "w") as f:
    f.write("10 08 08 //the magic number for this format file\n")
    f.write("#\n")
    f.write("# Smooth resolving power for Compton RT code.\n")
    f.write("# R(E) = R_base + (R_peak - R_base) * exp(-0.5*z^2)\n")
    f.write("# where z = (log10(E) - log10(E_center)) / sigma_dex\n")
    f.write(f"# R_base={R_base:g}  R_peak={R_peak:g}  "
            f"E_center={E_center:g} keV  sigma={sigma:g} dex\n")
    f.write(f"# {N_zones} zones from {E_lo_eV:g} eV to {E_hi_keV:g} keV, "
            f"smooth transitions\n")
    f.write("#\n")
    f.write("# Format: upper_limit_Ryd  resolving_power\n")
    f.write("# Last entry: 0 = upper bound of code\n")
    f.write("#\n")
    for e, r in zip(E_Ryd, R):
        f.write(f"{e:<14.4f} {r:.1f}\n")
    f.write(f"0        {closing_R:.6f}\n")

print(f"Wrote smooth_hump.ini with {len(E_Ryd)} entries")

# --- Visualize the resolution ---
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))

ax[0].plot(E_keV, R, "o-", ms=3, lw=1.2, color="C0")
ax[0].axvline(E_center, ls="--", color="gray", alpha=0.6,
              label=f"E_center = {E_center:g} keV")
ax[0].axhline(R_base, ls=":", color="gray", alpha=0.6,
              label=f"R_base = {R_base:g}")
ax[0].axhline(R_peak, ls=":", color="C3", alpha=0.6,
              label=f"R_peak = {R_peak:g}")
ax[0].set_xscale("log")
ax[0].set_xlabel("Energy [keV]")
ax[0].set_ylabel("Resolving power R = E/ΔE")
ax[0].set_title("Smooth-hump resolution")
ax[0].grid(True, which="both", alpha=0.3)
ax[0].legend(fontsize=8)

ax[1].plot(E_keV, 1.0 / R, "o-", ms=3, lw=1.2, color="C3")
ax[1].set_xscale("log")
ax[1].set_yscale("log")
ax[1].set_xlabel("Energy [keV]")
ax[1].set_ylabel("ΔE / E = 1/R")
ax[1].set_title("Relative bin width")
ax[1].grid(True, which="both", alpha=0.3)

plt.tight_layout()
plt.show()
plt.savefig("smooth_hump_resolution.png", dpi=150)
print("Saved smooth_hump_resolution.png")
