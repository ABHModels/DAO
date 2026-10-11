#ifndef DAO_NTDISK_H
#define DAO_NTDISK_H

#include "constants.h"

// Standalone, zero-ISCO-stress Novikov--Thorne annuli. No RT/Cloudy state.
// Flux and temperature are local comoving quantities, not observer quantities.
namespace ntdisk {

// Fixed diagnostic composition: unscaled Cloudy solar84.abn. Mass per H is
// sum_Z abundance(Z)*AtomicWeight(Z)*amu, using the Cloudy dense.h weights and
// amu=1.66053906660e-24 g (30 elements, H through Zn). No runtime Cloudy call.
constexpr double reference_mass_per_hydrogen_g = 2.3792220631064633e-24;
// User's reference convention n_e=1.21*n_H, rho=mass_per_H*n_H:
// kappa_T=sigma_T*n_e/rho. This is a composition conversion, not HH98's pure-H
// sigma_e/m_H. The electron ratio matches DAO's reference depth coordinate.
constexpr double reference_thomson_opacity_cm2_g =
    phys::reference_electrons_per_hydrogen * phys::sigma_T / reference_mass_per_hydrogen_g;

struct SystemParameters {
    double mass_msun;               // Black-hole mass [solar masses].
    double accretion_rate_g_s;       // Physical rest-mass accretion rate [g/s].
    double spin;                    // a*=cJ/(GM^2), -1<a*<1; disk rotates positively.
};

struct AnnulusRadiation {
    double radius_rg;               // Boyer--Lindquist R/(GM/c^2).
    double isco_rg;
    double intrinsic_flux_erg_cm2_s; // One face; no irradiation/reprocessing.
    double effective_temperature_K;
};

struct DiagnosticParameters {
    double alpha0;                  // Vertically averaged stress / TOTAL pressure.
    // Default uses the solar84, n_e=1.21*n_H reference opacity above.
    // HH98's pure-H sigma_e/m_H is explicitly generalized to this constant.
    // An override changes the analytic Sigma, not the reference tau conversion.
    double electron_opacity_cm2_g = reference_thomson_opacity_cm2_g;
};

struct SurfaceDensityDiagnostic {
    // Radiation-pressure / electron-scattering approximation, not a full disk
    // structure solution or a check that a chosen annulus is in that regime.
    double surface_density_g_cm2;   // BOTH faces: Sigma=2*m0.
    double mean_pressure_over_density_cm2_s2;
    double thomson_depth_midplane;  // Surface to midplane: kappa_T*Sigma/2.
    double thomson_depth_total;     // Across both faces: kappa_T*Sigma.
};

// Bardeen, Press & Teukolsky (1972), ApJ 178, 347, Eq. (2.21).
double isco_radius_rg(double spin);

// Hubeny & Hubeny (1998), ApJ 505, 558, astro-ph/9804288, Eq. (10).
// R>=ISCO; exactly at ISCO the zero-stress flux and temperature are zero.
// alpha0 and Sigma cannot affect this interface: neither is an input.
AnnulusRadiation evaluate_radiation(const SystemParameters&, double radius_rg);

// Hubeny & Hubeny (1998), Eqs. (14)--(18); R>ISCO required.
// Diagnostic only. Never supplies density, column, height or pressure to RT.
SurfaceDensityDiagnostic evaluate_surface_density_diagnostic(
    const SystemParameters&, double radius_rg, const DiagnosticParameters&);

} // namespace ntdisk
#endif
