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
#include "hmi.h"
#include "lines.h"
#include "transition.h"
#include "iso.h"
#include "taulines.h"
#include "physconst.h"

void CloudyInput::init(const ModelParams& par)
{
	nh  = par.nh;
	Afe = par.Afe;
}

// Pull emissivity and opacity from Cloudy after cdDrive(). Cloudy runs as a
// single-zone atomic-physics backend ("no line transfer", "stop zone 1");
// our RT solver does the transfer. Emissivity has three components:
// (A) continuum, (B) bound-bound lines, (C) inner-shell fluorescence.
// Opacity (kabs) is continuum only -- line self-absorption is handled by the
// per-line escape probability in (B).
void extract_cloudy_output(int id, RadField& rad, const RTGrids& g, int outer_iter)
{
	// Thermodynamic state
	rad.T_K[id]     = cdTemp_last();
	rad.n_e[id]     = cdEDEN_last();
	rad.heating[id] = cdHeating_last();
	rad.cooling[id] = cdCooling_last();

	// Map our energy grid to Cloudy's continuum index (j0 = Cloudy index of our bin 0)
	double E_lo_ryd = g.ene[0] / phys::eV_per_Ryd;
	long j0 = 0;
	while (j0 < rfield.nflux && rfield.anu(j0) < E_lo_ryd * 0.999)
		++j0;

	for (int i = 0; i < g.NE; ++i)
	{
		rad.jnu[id][i]      = 0.0;
		rad.kabs[id][i]     = 0.0;
		rad.jnu_line[id][i] = 0.0;
	}

	// (A) Continuum emissivity + opacity. ConEmitLocal is photons cm^-3 s^-1
	// per bin; convert to erg cm^-3 s^-1 eV^-1 sr^-1 via anu*eV_to_erg/widflx/4pi.
	// opacity_abs = continuum absorption (photoionization, free-free, H-, ...).
	for (int i = 0; i < g.NE; ++i)
	{
		long j = j0 + i;
		double conv1 = rfield.anu(j) * phys::eV_to_erg
		             / rfield.widflx(j) / phys::four_pi;

		rad.jnu[id][i]  = double(rfield.ConEmitLocal[nzone][j]) * conv1;
		rad.kabs[id][i] = opac.opacity_abs[j];
	}

	// Remove Cloudy's bound-electron Compton recoil opacity from kabs.
	// Cloudy folds it into opacity_abs as absorption (opacity_addtotal.cpp), and
	// "No scattering opacity" does not remove it. Our solver handles all electron
	// scattering via ksct + the Compton kernel, so leaving it would double-count
	// scattering as absorption (crashing the 20-100 keV albedo and erasing the
	// Compton hump). Subtract exactly what Cloudy added so kabs is pure photoabs.
	for (long nelem = 0; nelem < LIMELM; ++nelem)
	{
		if (!dense.lgElmtOn[nelem]) continue;
		for (long ion = 0; ion < nelem + 1; ++ion)
		{
			double factor = dense.xIonDense[nelem][ion];
			if (nelem == ipHYDROGEN) factor += hmi.H2_total * 2.0;
			if (factor <= 0.0) continue;
			factor *= ionbal.nCompRecoilElec[nelem - ion];
			if (factor <= 0.0) continue;
			long jstart = ionbal.ipCompRecoil[nelem][ion] - 1;
			for (int i = 0; i < g.NE; ++i)
			{
				long j = j0 + i;
				if (j < jstart || j - 1 + opac.iopcom < 0) continue;
				rad.kabs[id][i] -= opac.OpacStack[j - 1 + opac.iopcom] * factor;
			}
		}
	}
	for (int i = 0; i < g.NE; ++i)
		if (rad.kabs[id][i] < 0.0) rad.kabs[id][i] = 0.0;

	// (B) Bound-bound lines with per-line escape probability beta(tau_cell),
	// where beta depends on the redistribution type (PRD / CRD-wing / CRD-core).
	// "no line transfer" leaves Pesc=1 for all lines except Ly-alpha, which
	// Cloudy still computes; undo that Pesc first, then apply our own beta.
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

		// Undo Cloudy's residual Ly-alpha Pesc
		realnum Ptot = tr.Emis().Pesc_total();
		if (tr.Emis().iRedisFun() == ipLY_A && Ptot > 0. && Ptot < 1.)
			emiss /= Ptot;

		// Per-line escape probability beta(tau_cell)
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

		// Escaped line emission -> jnu_line (sets spectral shape)
		int i = (int)(j - j0);
		double dE_eV = rfield.widflx(j) * phys::eV_per_Ryd;
		double emiss_per_eV_sr = emiss / (dE_eV * phys::four_pi);
		rad.jnu_line[id][i] += emiss_per_eV_sr*beta;
	}

	// (C) Inner-shell fluorescence (e.g. Fe K-alpha 6.4 keV), no escape probability.
	// A K-shell photoionization leaves a vacancy that decays radiatively with
	// fluorescence yield omega (Auger otherwise). The photon lies below the
	// daughter ion's K-edge, so it is not reabsorbed by that transition -- the
	// only sink is continuum K-shell photoionization (already in opacity_abs).
	// Rate n_phot = n_ion * Gamma_shell * omega, following Cloudy's prt_lines.cpp;
	// atomic data from t_yield (Kaastra & Mewe 1993, A&AS 97, 443).
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

	// Combine: jnu = continuum (A) + escaped lines (B) + fluorescence (C)
	for (int i = 0; i < g.NE; ++i)
		rad.jnu[id][i] += rad.jnu_line[id][i];
}

// Cloudy commands issued once before the depth loop. Cloudy is run as a
// single-zone atomic-physics backend with its own RT disabled; our solver
// handles continuum attenuation, Compton scattering, and line escape.
void CloudyInput::issue_constant()
{
	cdRead("abundances \"solar84.abn\"");  // Grevesse & Anders 1984 solar
	cdRead("init \"xray12.ini\"");         // X-ray optimised continuum mesh
	char abuf[256];
	snprintf(abuf, sizeof(abuf), "element scale iron %.4f", Afe);
	cdRead(abuf);
	snprintf(abuf, sizeof(abuf), "hden %.6f", nh);
	cdRead(abuf);

	// Disable physics not needed for X-ray RT
	cdRead("no grain physics");
	cdRead("no fine opacities");
	cdRead("no molecules");
	cdRead("No scattering opacity");   // we handle Compton scattering ourselves
	cdRead("no level2 lines");

	// "no line transfer" sets Pesc=1 for all lines except Ly-alpha, which we
	// undo in extract_cloudy_output; we apply our own per-line beta(tau).
	cdRead("no line transfer");
	cdRead("Database H-like Lyman pumping off");

	cdRead("iterate convergence");
	cdRead("High temperature approach");   // numerical stability at T > 1e7 K
	cdRead("stop zone 1");
	cdRead("set temperature convergence 0.005");

	cdRead("Database Chianti mixed");   // CHIANTI where available, else Stout
}

// Write the incident SED as a Cloudy "table SED" and set its intensity.
// SED column is J0 directly (= F_nu shape); do NOT apply wid/E (= log-grid
// spacing, which would distort the shape). "units eV" must be on the first
// data line. Values are floored at 1e-30 (not dropped) so absorption troughs
// are preserved rather than bridged over by Cloudy's log-log interpolation.
void CloudyInput::issue_depth(int id, const RadField& rad, const RTGrids& g,const ModelParams& par)
{
	FILE* fsed = open_data("SED_TEST_API_INCI.dat", "w");
	fprintf(fsed, "# E_eV  J0 (F_nu shape)\n");

	bool units_written = false;
	for (int i = 0; i < g.NE; ++i)
	{
		double val = rad.J0[id][i];
		if (val < 1e-30) val = 1e-30;
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

	char buf[256];
	snprintf(buf, sizeof(buf), "intensity %.8f range %.8f to %.8f ev",(rad.log_xi[id]-log10(phys::four_pi))+par.nh,g.E_IN_LO,g.E_IN_HI);
	cdRead(buf);
}

// Bootstrap Cloudy once to extract its energy grid.
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

// Final post-convergence pass: same SED as issue_depth(), plus dump iron fractions.
void CloudyInput::issue_depth_lastest(int id, const RadField& rad, const RTGrids& g,const ModelParams& par)
{
	FILE* fsed = open_data("SED_TEST_API_INCI.dat", "w");
	fprintf(fsed, "# E_eV  J0 (F_nu shape)\n");

	bool units_written = false;
	for (int i = 0; i < g.NE; ++i)
	{
		double val = rad.J0[id][i];
		if (val < 1e-30) val = 1e-30;
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

	char buf[256];
	snprintf(buf, sizeof(buf), "intensity %.8f range %.8f to %.8f ev",(rad.log_xi[id]-log10(phys::four_pi))+par.nh,g.E_IN_LO,g.E_IN_HI);
	cdRead(buf);

	char fname[256];
	snprintf(fname, sizeof(fname), "\"%s%i.iron\"", par.run_hash,id);
	snprintf(buf, sizeof(buf), "save element iron %s", fname);
	cdRead(buf);
}