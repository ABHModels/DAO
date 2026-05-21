#include "cloudy_interface.h"
#include "constants.h"
#include "cddefines.h"
#include "cddrive.h"
#include "rfield.h"
#include "opacity.h"
#include "compton_cross_section.h"
#include "rt_escprob.h"
#include "doppvel.h"
#include <cmath>
#include <cstdio>
#include <vector>
#include <sys/stat.h>
#include "radius.h"
#include "yield.h"
#include "ionbal.h"
#include "dense.h"
#include "lines.h"
#include "transition.h"
#include "iso.h"
#include "taulines.h"
#include "physconst.h"

// ============================================================
void CloudyInput::init(const ModelParams& par)
{
	nh  = par.nh;
	Afe = par.Afe;
}

// ============================================================
// extract_cloudy_output
//
// Pull emissivity and opacity from Cloudy after cdDrive().
//
// Cloudy runs as a thin single-zone atomic physics backend
// ("no line transfer", "stop zone 1").  It provides level
// populations, ionization balance, and cross-sections.
// Our RT solver handles the actual radiative transfer.
//
// The emissivity has three components:
//
//   (A) Continuum (free-free, free-bound, two-photon)
//       Source: rfield.ConEmitLocal
//       No escape probability — smooth continuum, no self-absorption.
//
//   (B) Bound-bound line emission
//       Source: cdEmis_ip (Hazy2 Section 8.6.6)
//       Each line gets its own escape probability beta(tau):
//         emiss_eff = emiss_raw * beta(tau_cell)
//       where tau_cell = kappa_line * dr_cell.
//       This accounts for line self-absorption within the RT cell.
//       The escape function beta depends on the line's redistribution
//       type (PRD for resonance lines, CRD for others).
//
//   (C) Inner-shell fluorescence (e.g. Fe K-alpha 6.4 keV)
//       Source: t_yield (Kaastra & Mewe 1993)
//       No escape probability — the fluorescence photon energy is
//       below the K-edge of the daughter ion, so it cannot be
//       reabsorbed by the same bound-bound transition.  Absorption
//       is by K-shell photoionization of other ions (continuum
//       opacity, already in opacity_abs).
//
// The opacity is continuum only:
//   kabs = opacity_abs (photoionization + free-free + H- + etc.)
//   Line opacity is NOT included — line self-absorption is handled
//   by the per-line escape probability above.
//
// Note on Ly-alpha:
//   With "no line transfer", Cloudy still computes escape probability
//   for Ly-alpha lines (iRedisFun == ipLY_A) of all iso-sequence
//   species (rt_line_one.cpp:410).  So cdEmis_ip for Ly-alpha returns
//   emiss * Pesc_cloudy.  We divide out Pesc_cloudy before applying
//   our own beta(tau_cell) to avoid double-counting.
//   For all other lines, Pesc_cloudy = 1 (no correction needed).
//
void extract_cloudy_output(int id, RadField& rad, const RTGrids& g, int outer_iter)
{
	// ---- Thermodynamic state from Cloudy ----
	rad.T_K[id]     = cdTemp_last();     // electron temperature [K]
	rad.n_e[id]     = cdEDEN_last();     // electron density [cm^-3]
	rad.heating[id] = cdHeating_last();  // total heating [erg cm^-3 s^-1]
	rad.cooling[id] = cdCooling_last();  // total cooling [erg cm^-3 s^-1]

	// ---- Map our energy grid to Cloudy's continuum index ----
	// j0 is the Cloudy index corresponding to our bin 0.
	double E_lo_ryd = g.ene[0] / phys::eV_per_Ryd;
	long j0 = 0;
	while (j0 < rfield.nflux && rfield.anu(j0) < E_lo_ryd * 0.999)
		++j0;

	// ---- Zero output arrays ----
	for (int i = 0; i < g.NE; ++i)
	{
		rad.jnu[id][i]      = 0.0;
		rad.kabs[id][i]     = 0.0;
		rad.jnu_line[id][i] = 0.0;
	}

	// =============================================================
	// (A) Continuum emissivity and opacity
	//
	// ConEmitLocal[nzone][j]: photons cm^-3 s^-1 per energy bin.
	// Convert to erg cm^-3 s^-1 eV^-1 sr^-1:
	//   multiply by h*nu [erg/photon] = anu [Ryd] * eV_to_erg
	//   divide by bin width [eV]      = widflx [Ryd] * eV_per_Ryd
	//   divide by 4*pi [sr]
	// The eV_per_Ryd cancels, giving: anu * eV_to_erg / widflx / 4pi
	//
	// opacity_abs: continuum absorption [cm^-1]
	//   Includes photoionization, free-free, H-, H2+, grains.
	//   Does NOT include line opacity.
	// =============================================================
	for (int i = 0; i < g.NE; ++i)
	{
		long j = j0 + i;
		double conv1 = rfield.anu(j) * phys::eV_to_erg
		             / rfield.widflx(j) / phys::four_pi;

		rad.jnu[id][i]  = double(rfield.ConEmitLocal[nzone][j]) * conv1;
		rad.kabs[id][i] = opac.opacity_abs[j];
	}

	// =============================================================
	// (B) Bound-bound line emission with per-line escape probability
	//
	// For each line:
	//   1. Get raw emissivity from cdEmis_ip
	//   2. Undo Cloudy's Pesc for Ly-alpha lines
	//   3. Compute per-line κ, τ, β(τ)
	//   4. Escaped emission: emiss × β  → goes into jnu_line (spectral)
	// =============================================================
	double dr_cell = g.dr[id];

	for (long ip = 0; ip < LineSave.nsum; ++ip)
	{
		LinSv& line = LineSave.lines[ip];

		if (line.chSumTyp() == 'i') continue;
		if (line.isBlend()) continue;

		double emiss;
		cdEmis_ip(ip, &emiss, false);
		if (emiss <= 0.) continue;

		TransitionProxy tr = line.getTransition();
		if (!tr.associated()) continue;
		if (tr.ipCont() <= 0) continue;
		long j = tr.ipCont() - 1;
		if (j < j0 || j >= j0 + g.NE) continue;

		// Undo Cloudy's Pesc for Ly-alpha type lines
		realnum Ptot = tr.Emis().Pesc_total();
		if (tr.Emis().iRedisFun() == ipLY_A && Ptot > 0. && Ptot < 1.) 
			emiss /= Ptot;

		// Per-line escape probability
		double beta = 1.0;
		double pop = tr.Emis().PopOpc();
		double op  = tr.Emis().opacity();
		if (pop > 0. && op > 0.)
		{
			double kappa = pop * op * rfield.anu(j)
			             / (SPEEDLIGHT * rfield.widflx(j));
			if (kappa > 0.)
			{
				double tau_cell = kappa * dr_cell;
				int redis = tr.Emis().iRedisFun();
				if (redis == ipPRD || redis == ipLY_A)
				{
					double damp = tr.Emis().dampXvel()
					            / GetDopplerWidth(dense.AtomicWeight[(*tr.Hi()).nelem()-1]);
					beta = esc_PRD_1side(tau_cell, damp);
				}
				else if (redis == ipCRDW)
				{
					double damp = tr.Emis().dampXvel()
					            / GetDopplerWidth(dense.AtomicWeight[(*tr.Hi()).nelem()-1]);
					beta = esc_CRDwing_1side(tau_cell, damp);
				}
				else if (redis == ipCRD)
				{
					beta = esca0k2(tau_cell);
				}

			}
		}

		// Escaped line emission → jnu_line (determines spectral shape)
		int i = (int)(j - j0);
		double dE_eV = rfield.widflx(j) * phys::eV_per_Ryd;
		double emiss_per_eV_sr = emiss / (dE_eV * phys::four_pi);
		rad.jnu_line[id][i] += emiss_per_eV_sr*beta;
	}

	// =============================================================
	// (C) Inner-shell fluorescence (no escape probability)
	//
	// Physical mechanism:
	//   1. An X-ray photon (E > K-edge) ejects a K-shell electron:
	//        X^+q + gamma -> X^+(q+1) + e^-(K-shell)
	//
	//   2. An outer-shell electron fills the K-shell vacancy.
	//      Two competing decay channels:
	//        Auger:         energy ejects another electron (no photon)
	//        Fluorescence:  energy emitted as an X-ray photon
	//      The fluorescence yield omega = P(photon) / P(total).
	//      omega increases with Z:  ~0.01 (O), ~0.34 (Fe), ~0.97 (U)
	//      because the radiative rate scales as Z^4 while the Auger
	//      rate is roughly constant.
	//
	//   3. The fluorescence photon (e.g. Fe K-alpha at 6.4 keV)
	//      has energy BELOW the K-edge of the daughter ion.
	//      It cannot be reabsorbed by the same bound-bound transition
	//      — there is no such transition.  The only absorption channel
	//      is K-shell photoionization of less-ionized ions, which is
	//      continuum opacity (already in opacity_abs).
	//      Therefore: no escape probability is needed.
	//
	// Data source: Kaastra & Mewe (1993, A&AS 97, 443).
	// Cloudy stores the atomic data in t_yield.
	//
	// Rate:  n_phot = n_ion * Gamma_shell * omega  [phot cm^-3 s^-1]
	//   n_ion       = parent ion density [cm^-3]
	//   Gamma_shell = inner-shell photoionization rate [s^-1]
	//   omega       = fluorescence yield (dimensionless, 0 to 1)
	// =============================================================
	for (long ifl = 0; ifl < t_yield::Inst().nlines(); ++ifl)
	{
		long ip = t_yield::Inst().ipoint(ifl) - 1;
		if (ip < j0 || ip >= j0 + g.NE) continue;

		double n_phot =
			dense.xIonDense[t_yield::Inst().nelem(ifl)][t_yield::Inst().ion(ifl)] *
			ionbal.PhotoRate_Shell[t_yield::Inst().nelem(ifl)]
								  [t_yield::Inst().ion(ifl)]
								  [t_yield::Inst().nshell(ifl)][0] *
			t_yield::Inst().yield(ifl);
		if (n_phot <= 0.) continue;

		int i = (int)(ip - j0);
		double conv_fl = rfield.anu(ip) * phys::eV_to_erg / rfield.widflx(ip);
		rad.jnu_line[id][i] += n_phot * conv_fl / phys::four_pi;
	}

	// =============================================================
	// Combine: jnu_total = jnu_continuum (boosted) + jnu_line (escaped)
	//
	// Per-line Pesc was already applied in (B) above.
	// Fluorescence in (C) has no Pesc (physically correct).
	// Continuum boost in (D) conserves energy.
	// =============================================================
	for (int i = 0; i < g.NE; ++i)
		rad.jnu[id][i] += rad.jnu_line[id][i];
}

// ============================================================
// issue_constant — Cloudy commands that stay the same for all
// depth points.  Called once before the depth loop.
//
// Design philosophy:
//   Cloudy runs as a thin single-zone ATOMIC PHYSICS backend.
//   It provides level populations, ionization balance, emissivity,
//   and opacity.  Our external RT solver handles all radiative
//   transfer (continuum attenuation, Compton scattering, line
//   escape probability).  Therefore we disable Cloudy's own RT
//   and simplify its physics to what we need.
//
void CloudyInput::issue_constant()
{
	// --- Abundances and atomic data ---
	cdRead("abundances \"solar84.abn\"");  // Grevesse & Anders 1984 solar
	cdRead("init \"xray12.ini\"");         // X-ray optimised continuum mesh
	char abuf[256];
	snprintf(abuf, sizeof(abuf), "element scale iron %.4f", Afe);
	cdRead(abuf);                          // scale Fe abundance
	snprintf(abuf, sizeof(abuf), "hden %.6f", nh);
	cdRead(abuf);                          // hydrogen density [cm^-3]

	// --- Disable physics not needed for X-ray RT ---
	cdRead("no grain physics");       // no dust (X-ray regime)
	cdRead("no fine opacities");      // no fine-resolution continuum opacity
	cdRead("no molecules");           // no H2, CO, etc. (hot plasma)
	cdRead("No scattering opacity");  // we handle Compton scattering ourselves
	cdRead("no level2 lines");        // disable minor lines (speed)
	
	// --- Line transfer: handled by OUR RT solver ---
	// "no line transfer" makes Cloudy set Pesc=1 for all lines
	// EXCEPT Ly-alpha (iRedisFun == ipLY_A), which Cloudy still
	// computes internally (rt_line_one.cpp:410).
	// We apply our own per-line escape probability beta(tau_cell)
	// in extract_cloudy_output, and undo Cloudy's Pesc for Ly-alpha
	// to avoid double-counting.
	cdRead("no line transfer");

	// --- Disable H-like Ly-alpha pumping from the incident SED ---
	cdRead("Database H-like Lyman pumping off");

	// --- Solver settings ---
	cdRead("iterate convergence");          // iterate until converged
	cdRead("High temperature approach");    // numerical stability at T > 10^7 K
	cdRead("stop zone 1");                  // single thin zone (atomic physics only)
	cdRead("set temperature convergence 0.005");  // 0.5% T convergence

	// --- Atomic database ---
	// "mixed" = use CHIANTI where available, Stout for the rest.
	// CHIANTI provides better data for highly ionized species
	// (Fe XVII–XXIV, Si XIII, S XV, etc.).
	cdRead("Database Chianti mixed");
}

// ============================================================
void CloudyInput::issue_depth(int id, const RadField& rad, const RTGrids& g,const ModelParams& par)
{
	// --- SED shape: write J0/E as table SED ---
	// Cloudy's "table SED" expects: E  SED(E)
	// where SED(E) = J0(E) × ΔE / E  (photon-number-weighted shape).
	// "units eV" MUST appear on the first data line — otherwise
	// Cloudy reads energies as Rydbergs (default), causing a
	// factor-of-13.6 energy mismatch.
	FILE* fsed = open_data("SED_TEST_API_INCI.dat", "w");
	fprintf(fsed, "# E_eV  J0*wid/E\n");

	bool units_written = false;
	for (int i = 0; i < g.NE; ++i)
	{
		double val = rad.J0[id][i] * g.wid[i] / g.ene[i];
		if (val <= 1e-30) continue;
		if (!units_written)
		{
			fprintf(fsed, "%.8e  %.8e units eV\n", g.ene[i], val);
			units_written = true;
		}
		else
		{
			fprintf(fsed, "%.8e  %.8e\n", g.ene[i], val);
		}
	}
	fprintf(fsed, "************");
	fclose(fsed);
	cdRead("table SED \"SED_TEST_API_INCI.dat\"");

	// --- Intensity normalisation ---
	char buf[256];
	snprintf(buf, sizeof(buf), "intensity %.8f range %.8f to %.8f ev",(rad.log_xi[id]-log10(phys::four_pi))+par.nh,g.E_IN_LO,g.E_IN_HI);
	cdRead(buf);
}

// ============================================================
// Bootstrap Cloudy to extract the energy grid.
// ============================================================
void bootstrap_cloudy_energy_grid(RTGrids& g, const char* save_file)
{
	fprintf(stdout, "Bootstrapping Cloudy for energy grid...\n");

	cdInit();
	cdTalk(false);
	cdRead("init \"xray12.ini\"");
	cdRead("hden 16");
	cdRead("AGN 6.0 -1.4 -0.5 -1.0");
	cdRead("intensity 5.0");
	cdRead("radius 10 -2");
	cdRead("no level2 lines");
	cdRead("no molecules");
	cdRead("no grain physics");
	cdRead("stop zone 1");
	cdDrive();

	int nf = rfield.nflux;
	std::vector<double> anu(nf), widf(nf);
	for (int j = 0; j < nf; ++j)
	{
		anu[j]  = rfield.anu(j);
		widf[j] = rfield.widflx(j);
	}
	g.init_energy_from_cloudy(nf, anu.data(), widf.data());
	g.save_energy(save_file);

	fprintf(stdout, "  Energy grid: NE=%d  E=[%.3f, %.3f] eV\n",
	        g.NE, g.ene[0], g.ene[g.NE - 1]);
			
	char buf[256];
	cdVersion(buf);
	fprintf(stdout, "  Cloudy version: %s\n",buf);
}

void CloudyInput::issue_depth_lastest(int id, const RadField& rad, const RTGrids& g,const ModelParams& par)
{
	// after the whole iteration, output the line labels and the iron fraction
	FILE* fsed = open_data("SED_TEST_API_INCI.dat", "w");
	fprintf(fsed, "# E_eV  J0*wid/E\n");

	bool units_written = false;
	for (int i = 0; i < g.NE; ++i)
	{
		double val = rad.J0[id][i] * g.wid[i] / g.ene[i];
		if (val <= 1e-30) continue;
		if (!units_written)
		{
			fprintf(fsed, "%.8e  %.8e units eV\n", g.ene[i], val);
			units_written = true;
		}
		else
		{
			fprintf(fsed, "%.8e  %.8e\n", g.ene[i], val);
		}
	}
	fprintf(fsed, "************");
	fclose(fsed);
	cdRead("table SED \"SED_TEST_API_INCI.dat\"");

	// --- Intensity normalisation ---
	char buf[256];
	snprintf(buf, sizeof(buf), "intensity %.8f range %.8f to %.8f ev",(rad.log_xi[id]-log10(phys::four_pi))+par.nh,g.E_IN_LO,g.E_IN_HI);
	cdRead(buf);

	// line labels
	char fname[256];
	snprintf(fname, sizeof(fname), "\"%s%i.iron\"", par.run_hash,id);
	snprintf(buf, sizeof(buf), "save element iron %s", fname);
	cdRead(buf);
}