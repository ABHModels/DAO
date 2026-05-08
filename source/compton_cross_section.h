#ifndef COMPTON_CROSS_SECTION_H
#define COMPTON_CROSS_SECTION_H

#include "constants.h"

// ============================================================
// Compton scattering cross section σ(E, T)
//
// Accounts for:
//   - Klein-Nishina corrections at high photon energies
//   - Relativistic Maxwellian electron distribution at high T
//
// Based on Poutanen & Svensson (1996)
// files on this repo: bk2.f,crsexact.f,scattxs.f
// which implements the exact cross section from
// Poutanen & Svensson (1996).
//
// Three regimes:
//   1. θ < 0.0169 (T < ~10^8 K) and x < 0.002 (E < ~1 keV):
//      σ = σ_T  (Thomson limit)
//   2. θ < 0.0169 and x ≥ 0.002:
//      σ = Klein-Nishina total cross section
//   3. θ ≥ 0.0169 (T ≥ ~10^8 K):
//      σ = exact relativistic, averaged over Maxwellian electrons
//
// where x = E / (m_e c^2),  θ = kT / (m_e c^2)
// ============================================================

// Single cross section [cm^2] for photon energy E_eV [eV]
// in a plasma at temperature T_K [K].
double compton_cross_section(double E_eV, double T_K);

// Fill scattering opacity array for all energies at one temperature:
//   ksct[ie] = n_e * σ(ene_eV[ie], T_K)
void compute_compton_opacity(double* ksct, int NE, const double* ene_eV,
                             double T_K, double n_e);

#endif // COMPTON_CROSS_SECTION_H
