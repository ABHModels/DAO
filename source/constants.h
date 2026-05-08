#ifndef CONSTANTS_H
#define CONSTANTS_H

// ============================================================
// Physical constants for the Compton scattering RT code
//
// CGS units unless noted otherwise.
// ============================================================

namespace phys {

// Fundamental constants
constexpr double h       = 6.62607015e-27;   // Planck constant [erg s]
constexpr double c       = 2.99792458e10;    // speed of light [cm s^-1]
constexpr double k_B     = 1.380649e-16;     // Boltzmann constant [erg K^-1]
constexpr double m_e     = 9.1093837015e-28; // electron mass [g]
constexpr double sigma_T = 6.6524587321e-25; // Thomson cross section [cm^2]

// Unit conversions
constexpr double eV_to_erg  = 1.602176634e-12;  // 1 eV in erg
constexpr double eV_per_Ryd = 13.6058;          // 1 Rydberg in eV
constexpr double keV_to_eV  = 1.0e3;            // 1 keV in eV
constexpr double m_e_c2_eV  = 5.10999e5;        // electron rest energy [eV]

// Numerical constants
constexpr double four_pi = 4.0 * 3.14159265358979323846;

} // namespace phys

// Temperature cache grid (shared by KernelCache and ScatteringCache)
const double T_CACHE_LO = 1e4;    // K
const double T_CACHE_HI = 1e9;    // K
const int    N_T_CACHE  = 50;     // number of temperature grid points

// Kernel cache directory: reads from KERNEL_DIR env var, defaults to "kernel"
#include <cstdlib>
inline const char* kernel_cache_dir()
{
	const char* dir = std::getenv("KERNEL_DIR");
	return dir ? dir : "kernel";
}

#endif // CONSTANTS_H
