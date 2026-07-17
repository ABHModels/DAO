#include "cloudy_interface.h"
#include "constants.h"
#include "cddefines.h"
#include "cddrive.h"
#include "rfield.h"
#include "opacity.h"
#include "compton_cross_section.h"
#include "rt_escprob.h"
#include "doppvel.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <cstdio>
#include <vector>
#include <string>
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
#include <sys/stat.h>

#undef fopen

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
void extract_cloudy_output(int id, RadField& rad, const RTGrids& g, int outer_iter,
                           std::vector<LineRec>& recs)
{
	// Thermodynamic state
	rad.T_K[id]     = cdTemp_last();
	rad.n_e[id]     = cdEDEN_last();
	rad.heating[id] = cdHeating_last();
	rad.line_heat[id] = 0.0;
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

	// (B) Bound-bound lines. Cloudy gives the local photon production rate and
	// atomic opacity; we postpone escape until the whole slab is known, because
	// a line photon can be trapped by absorbing cells far from where it was born.
	for (long ip = 0; ip < LineSave.nsum; ++ip)
	{
		// One physical or bookkeeping line in Cloudy's line stack.
		LinSv& line = LineSave.lines[ip];

		// Skip labels and blends: they do not carry independent line photons.
		if (line.chSumTyp() == 'i') continue;
		if (line.isBlend()) continue;

		// Atomic transition u->l: populations, A_ul, opacity, line profile.
		TransitionProxy tr = line.getTransition();
		if (!tr.associated()) continue;
		// Keep only lines whose photon energy lies on our RT energy grid.
		if (tr.ipCont() <= 0) continue;
		long j = tr.ipCont() - 1;
		if (j < j0 || j >= j0 + g.NE) continue;

		// Thin line power: 4pi*j_line = n_u*A_ul*h*nu. With Cloudy's line
		// transfer disabled this is already the Pesc=1 rate, except H I Lya.
		double emiss = 0.0;
		cdEmis_ip(ip, &emiss, false);
		if (emiss < 0.0) emiss = 0.0;

		// Line-center opacity
		// kappa_nu0 = (n_l - n_u g_l/g_u) * sigma_nu0
		// sigma_nu0 = H(a,0) * sigma_v / v_D
		double pop = tr.Emis().PopOpc();	// net absorbing population
		double op  = tr.Emis().opacity();   // line opacity factor from oscillator strength
		double vth = GetDopplerWidth(dense.AtomicWeight[(*tr.Hi()).nelem()-1]); // Doppler velocity width
		double a = 0.0, kappaL = 0.0;
		if (pop > 0. && op > 0. && vth > 0.)
		{
			a      = tr.Emis().dampXvel() / vth;  // line-center Voigt parameter
			kappaL = VoigtH0(a) * pop * op / vth; // line-center absorption coefficient
		}

		// Electron scattering can move a trapped line photon out of the narrow core.
		// Cloudy's electron-escape branching uses the lower-level opacity scale,
		// n_l * sigma_v / v_D, not PopOpc: stimulated emission reduces net line
		// extinction but does not remove lower-level ions as resonant scatterers.

		double nlo = (*tr.Lo()).Pop();  // lower level population
		double kappaLo = (nlo > 0. && op > 0. && vth > 0.)
		               ? nlo * op / vth : 0.0; 
		// Pure absorbers emit nothing locally but still block deeper photons.
		if (emiss <= 0.0 && kappaL <= 0.0) continue;

		// Collisional quenching parameter y = C_ul/A_ul. Negative collision
		// strengths are Cloudy bookkeeping/no-data values, not physical rates.
		double Aul = tr.Emis().Aul();     // spontaneous decay rate
		double gu  = (*tr.Hi()).g();      // upper-level statistical weight
		double Col = tr.Coll().col_str(); // collision strength
		double cdsqte = dense.cdsqte;     // electron collisional factor
		double y = 0.0;
		if (Aul > 0.0 && gu > 0.0 && Col > 0.0 && cdsqte > 0.0)
			y = Col * cdsqte / (gu * Aul);
		if (!std::isfinite(y) || y < 0.0)
			y = 0.0;

		// H I Lya is still treated specially inside Cloudy although "no line transfer". 
		// Convert it back
		realnum Ptot = tr.Emis().Pesc_total();
		if (tr.Emis().iRedisFun() == ipLY_A && Ptot > 0. && Ptot < 1.)
		{
			double Ploss = tr.Emis().Ploss();
			emiss *= (Ploss + y) / (Ptot * (1.0 + y));
		}

		// Per-cell record consumed by apply_line_escape()
		LineRec rec;
		rec.bin    = (int)(j - j0);                       // RT energy bin of h*nu_ul
		rec.ip     = ip;                                  // line identity across cells
		rec.emiss  = emiss;                               // thin emissivity [erg cm^-3 s^-1]
		rec.dE_eV  = rfield.widflx(j) * phys::eV_per_Ryd; // continuum bin width [eV]
		rec.E_eV   = rfield.anu(j) * phys::eV_per_Ryd;    // photon energy [eV]
		rec.wave   = line.wavlVac();                      // vacuum wavelength [Angstrom]
		rec.kappaL = kappaL;                              // line-center opacity [cm^-1]
		rec.kappaLo = kappaLo;                            // lower-level-only opacity for Pelec
		rec.damp   = a;                                   // Voigt damping parameter
		rec.y      = y;                                   // C_ul/A_ul
		rec.redis  = tr.Emis().iRedisFun();               // PRD/CRD/CRD-wing/Lya
		rec.label  = line.label();                        // Cloudy label + wavelength
		rec.comment = line.chComment();                   // spectroscopic identification
		recs.push_back(rec);
	}

	// (C) Inner-shell fluorescence. A K-shell photoionization creates a vacancy;
	// fluorescence yield omega gives the radiative branch instead of Auger loss.
	// No line escape probability is applied: the emitted photon is below the
	// daughter K edge, so continuum opacity handles any later absorption.
	// Photon rate: n_phot = n_ion * Gamma_shell * omega.
	for (long ifl = 0; ifl < t_yield::Inst().nlines(); ++ifl)
	{
		long ip = t_yield::Inst().ipoint(ifl) - 1;
		if (ip < j0 || ip >= j0 + g.NE) continue;

		double n_phot =
			dense.xIonDense[t_yield::Inst().nelem(ifl)][t_yield::Inst().ion(ifl)] * // ion density
			ionbal.PhotoRate_Shell[t_yield::Inst().nelem(ifl)]                      // shell photoionization rate
								  [t_yield::Inst().ion(ifl)]
								  [t_yield::Inst().nshell(ifl)][0] *
			t_yield::Inst().yield(ifl);                                             // fluorescence yield
		if (n_phot <= 0.) continue;

		int i = (int)(ip - j0);
		double conv_fl = rfield.anu(ip) * phys::eV_to_erg / rfield.widflx(ip); // photon rate -> per-eV power
		rad.jnu_line[id][i] += n_phot * conv_fl / phys::four_pi;
	}

	// NOTE: jnu now holds continuum (A) only; jnu_line holds fluorescence (C).
	// Bound-bound escape (B) and the final fold jnu += jnu_line happen in
	// apply_line_escape() after the whole depth column has been collected.
}

// Probability that a line photon escapes to one boundary without another line
// absorption. The redistribution type controls how easily photons reach the
// wings, where the line opacity is lower.
static double esc_oneside(int redis, double tau, double damp)
{
	// No line optical depth to this boundary: free escape.
	if (tau <= 0.0) return 1.0;
	// PRD resonance lines keep memory of wing photons, so wing escape is easier.
	if (redis == ipPRD || redis == ipLY_A)
		return (damp > 0.0) ? esc_PRD_1side(tau, damp) : esca0k2(tau);
	// CRD-wing lines escape through rare re-emission into the Lorentz wings.
	if (redis == ipCRDW)
		return (damp > 0.0) ? esc_CRDwing_1side(tau, damp) : esca0k2(tau);
	// Pure Doppler-core CRD uses Hummer's K2 escape probability.
	if (redis == ipCRD)
		return esca0k2(tau);
	// Unknown redistribution: do not suppress the line silently.
	return 1.0;
}

static std::string tab_safe(std::string s)
{
	for (char& c : s)
		if (c == '\t' || c == '\n' || c == '\r')
			c = ' ';
	return s;
}

struct LineEscapeDiag
{
	long ip = -1;
	double E_eV = 0.0;
	double wave = 0.0;
	double thin = 0.0;
	double escaped = 0.0;
	double heat_dest = 0.0;
	double beta_w = 0.0;
	double pelec_w = 0.0;
	double pdest_w = 0.0;
	double y_w = 0.0;
	double tau_mid_w = 0.0;
	double tau_col_max = 0.0;
	std::string label;
	std::string comment;
};

static void write_line_escape_diagnostics(
	const RadField& rad,
	const RTGrids& g,
	const ModelParams& par,
	int iter,
	const std::vector<LineEscapeDiag>& diag,
	const std::vector<std::vector<double>>& thin_depth,
	const std::vector<std::vector<double>>& esc_depth,
	const std::vector<std::vector<double>>& heat_depth,
	const std::vector<std::vector<double>>& beta_depth_w,
	const std::vector<std::vector<double>>& pelec_depth_w,
	const std::vector<std::vector<double>>& pdest_depth_w,
	const std::vector<std::vector<double>>& y_depth_w)
{
	const char* dir = par.run_dir[0] ? par.run_dir : "results";
	mkdir("results", 0755);
	if (par.run_dir[0])
		mkdir(dir, 0755);

	std::vector<int> active;
	active.reserve(diag.size());
	for (int L = 0; L < (int)diag.size(); ++L)
		if (diag[L].thin > 0.0)
			active.push_back(L);

	std::sort(active.begin(), active.end(), [&](int a, int b)
	{
		return diag[a].escaped > diag[b].escaped;
	});

	char fname[512];
	snprintf(fname, sizeof(fname), "%s/line_escape_lines_iter%03d.dat", dir, iter);
	FILE* fp = fopen(fname, "w");
	if (fp)
	{
		fprintf(fp, "# Per-line escape summary  iter=%d\n", iter);
		fprintf(fp, "# Sorted by slab-integrated escaped line power.\n");
		fprintf(fp, "# Col 1: rank\n");
		fprintf(fp, "# Col 2: line_ip [Cloudy LineSave index]\n");
		fprintf(fp, "# Col 3: E_eV\n");
		fprintf(fp, "# Col 4: Cloudy wavelength label value\n");
		fprintf(fp, "# Col 5: thin_line_power [erg cm^-3 s^-1, depth-summed]\n");
		fprintf(fp, "# Col 6: escaped_line_power [erg cm^-3 s^-1, depth-summed]\n");
		fprintf(fp, "# Col 7: continuum_destroyed_heat [erg cm^-3 s^-1, depth-summed]\n");
		fprintf(fp, "# Col 8: escaped/thin\n");
		fprintf(fp, "# Col 9: heat_dest/thin\n");
		fprintf(fp, "# Col 10: mean_beta [thin-power weighted]\n");
		fprintf(fp, "# Col 11: mean_Pelec [thin-power weighted]\n");
		fprintf(fp, "# Col 12: mean_Pdest [thin-power weighted]\n");
		fprintf(fp, "# Col 13: mean_y [thin-power weighted]\n");
		fprintf(fp, "# Col 14: mean_tauT_emiss [thin-power weighted]\n");
		fprintf(fp, "# Col 15: max_column_tau_line\n");
		fprintf(fp, "# Col 16: Cloudy_label\n");
		fprintf(fp, "# Col 17: Cloudy_comment\n");
		for (int rank = 0; rank < (int)active.size(); ++rank)
		{
			const int L = active[rank];
			const LineEscapeDiag& d = diag[L];
			const double inv = (d.thin > 0.0) ? 1.0 / d.thin : 0.0;
			fprintf(fp,
			        "%d\t%ld\t%.9e\t%.9e\t%.9e\t%.9e\t%.9e\t"
			        "%.9e\t%.9e\t%.9e\t%.9e\t%.9e\t%.9e\t%.9e\t%.9e\t%s\t%s\n",
			        rank + 1, d.ip, d.E_eV, d.wave, d.thin, d.escaped,
			        d.heat_dest, d.escaped * inv, d.heat_dest * inv,
			        d.beta_w * inv, d.pelec_w * inv, d.pdest_w * inv,
			        d.y_w * inv, d.tau_mid_w * inv, d.tau_col_max,
			        tab_safe(d.label).c_str(), tab_safe(d.comment).c_str());
		}
		fclose(fp);
		fprintf(stdout, "  Saved: %s\n", fname);
	}

	std::vector<int> selected;
	std::vector<std::string> selected_reason;
	auto already_selected = [&](int L) -> bool
	{
		for (int old : selected)
			if (old == L) return true;
		return false;
	};
	auto add_unique = [&](int L, const char* reason)
	{
		if (L < 0 || already_selected(L)) return;
		selected.push_back(L);
		selected_reason.push_back(reason);
	};

	double max_thin = 0.0;
	for (int L : active)
		if (diag[L].thin > max_thin)
			max_thin = diag[L].thin;

	auto escaped_fraction = [&](int L) -> double
	{
		return (diag[L].thin > 0.0) ? diag[L].escaped / diag[L].thin : 0.0;
	};

	int strongest_escaped = -1;
	for (int L : active)
	{
		if (strongest_escaped < 0 ||
		    diag[L].escaped > diag[strongest_escaped].escaped)
			strongest_escaped = L;
	}
	add_unique(strongest_escaped, "strongest escaped line");

	int strongest_destroyed = -1;
	for (int L : active)
	{
		// Require real continuum destruction: without this floor a
		// continuum-transparent slab would silently label a heat_dest==0
		// line as the "strongest continuum-destroyed line". Dropping the
		// slot lets the unthresholded fallback below still fill four lines.
		if (already_selected(L) || diag[L].thin < 0.01 * max_thin ||
		    diag[L].heat_dest <= 0.0)
			continue;
		if (strongest_destroyed < 0 ||
		    diag[L].heat_dest > diag[strongest_destroyed].heat_dest)
			strongest_destroyed = L;
	}
	add_unique(strongest_destroyed, "strongest continuum-destroyed line");

	int strongly_suppressed = -1;
	double suppressed_score = -1.0;
	const double ionizing_line_floor_eV = 13.6;
	for (int L : active)
	{
		if (already_selected(L) || diag[L].thin < 0.02 * max_thin)
			continue;
		if (diag[L].E_eV < ionizing_line_floor_eV)
			continue;
		const double frac = escaped_fraction(L);
		if (frac >= 0.35)
			continue;
		const double score = diag[L].thin * (1.0 - frac);
		if (score > suppressed_score)
		{
			suppressed_score = score;
			strongly_suppressed = L;
		}
	}
	if (strongly_suppressed < 0)
	{
		for (int L : active)
		{
			if (already_selected(L) || diag[L].thin < 0.02 * max_thin)
				continue;
			const double frac = escaped_fraction(L);
			if (frac >= 0.35)
				continue;
			const double score = diag[L].thin * (1.0 - frac);
			if (score > suppressed_score)
			{
				suppressed_score = score;
				strongly_suppressed = L;
			}
		}
	}
	add_unique(strongly_suppressed, "strong line with strong escape attenuation");

	int nearly_surviving = -1;
	for (int L : active)
	{
		if (already_selected(L) || diag[L].thin < 0.02 * max_thin)
			continue;
		if (escaped_fraction(L) < 0.95)
			continue;
		if (nearly_surviving < 0 || diag[L].escaped > diag[nearly_surviving].escaped)
			nearly_surviving = L;
	}
	add_unique(nearly_surviving, "strong line with weak escape attenuation");

	for (int L : active)
	{
		if ((int)selected.size() >= 4) break;
		add_unique(L, "fallback strongest escaped line");
	}

	snprintf(fname, sizeof(fname), "%s/line_escape_selected_iter%03d.dat", dir, iter);
	fp = fopen(fname, "w");
	if (fp)
	{
		fprintf(fp, "# Selected line escape depth profiles  iter=%d\n", iter);
		fprintf(fp, "# Selection: representative escape-physics cases from the per-line summary.\n");
		for (int s = 0; s < (int)selected.size(); ++s)
		{
			const int L = selected[s];
			const LineEscapeDiag& d = diag[L];
			const double inv = (d.thin > 0.0) ? 1.0 / d.thin : 0.0;
			fprintf(fp,
			        "# Selected %d: %s; ip=%ld; escaped/thin=%.6e; "
			        "heat_dest/thin=%.6e; label=%s\n",
			        s + 1, selected_reason[s].c_str(), d.ip,
			        d.escaped * inv, d.heat_dest * inv,
			        tab_safe(d.label).c_str());
		}
		fprintf(fp, "# Col 1: selected_index\n");
		fprintf(fp, "# Col 2: depth index\n");
		fprintf(fp, "# Col 3: tau_mid [Thomson]\n");
		fprintf(fp, "# Col 4: line_ip [Cloudy LineSave index]\n");
		fprintf(fp, "# Col 5: E_eV\n");
		fprintf(fp, "# Col 6: thin_line_power [erg cm^-3 s^-1]\n");
		fprintf(fp, "# Col 7: escaped_line_power [erg cm^-3 s^-1]\n");
		fprintf(fp, "# Col 8: escaped/thin\n");
		fprintf(fp, "# Col 9: continuum_destroyed_heat [erg cm^-3 s^-1]\n");
		fprintf(fp, "# Col 10: mean_beta [thin-power weighted in this depth]\n");
		fprintf(fp, "# Col 11: mean_Pelec [thin-power weighted in this depth]\n");
		fprintf(fp, "# Col 12: mean_Pdest [thin-power weighted in this depth]\n");
		fprintf(fp, "# Col 13: mean_y [thin-power weighted in this depth]\n");
		fprintf(fp, "# Col 14: Cloudy_label\n");
		fprintf(fp, "# Col 15: Cloudy_comment\n");
		for (int s = 0; s < (int)selected.size(); ++s)
		{
			const int L = selected[s];
			for (int id = 0; id < g.ND_MID; ++id)
			{
				const double thin = thin_depth[L][id];
				const double inv = (thin > 0.0) ? 1.0 / thin : 0.0;
				fprintf(fp,
				        "%d\t%d\t%.9e\t%ld\t%.9e\t%.9e\t%.9e\t%.9e\t%.9e\t"
				        "%.9e\t%.9e\t%.9e\t%.9e\t%s\t%s\n",
				        s + 1, id, g.tau_mid[id], diag[L].ip, diag[L].E_eV,
				        thin, esc_depth[L][id], esc_depth[L][id] * inv,
				        heat_depth[L][id], beta_depth_w[L][id] * inv,
				        pelec_depth_w[L][id] * inv, pdest_depth_w[L][id] * inv,
				        y_depth_w[L][id] * inv, tab_safe(diag[L].label).c_str(),
				        tab_safe(diag[L].comment).c_str());
			}
		}
		fclose(fp);
		fprintf(stdout, "  Saved: %s\n", fname);
	}
}

// Second pass over the full slab: apply column-scale line escape.
// A bound-bound photon is repeatedly absorbed and re-emitted until one channel
// wins: line escape, electron-scattering escape, collisional quenching, or
// continuum absorption. Only the surviving fraction emiss*P is injected.
//
// tau_in/out = line optical depth from the emission cell to each boundary.
// beta       = 0.5 * [esc(tau_in) + esc(tau_out)]        line escape
// y          = C_ul / A_ul                               collisional quench
// Pdest      = continuum absorption during trapping
// Pelec      = Thomson scatter out of the line core
// H_Pdest   = emiss*(1+y)*Pdest/(beta + Pelec + y + Pdest)
// P          = (beta + Pelec)*(1+y)/(beta + Pelec + y + Pdest)
// j_line     = emiss * P / (dE * 4pi)
// H_Pdest is accumulated in rad.line_heat for the next Cloudy thermal pass.
void apply_line_escape(RadField& rad, const RTGrids& g,
                       const std::vector<std::vector<LineRec>>& store,
                       const ModelParams& par, int iter)
{
	const int ND = g.ND_MID;

	// The same transition appears in many depth cells; give each transition a
	// compact column index L.
	std::map<long,int> idx;
	for (int id = 0; id < ND; ++id)
		for (const LineRec& r : store[id])
			if (idx.find(r.ip) == idx.end())
			{
				int n = (int)idx.size();
				idx[r.ip] = n;
			}
	const int NL = (int)idx.size();
	// No lines inside the RT window: still fold fluorescence (C) into jnu
	if (NL == 0)
	{
		for (int id = 0; id < ND; ++id)
			for (int i = 0; i < g.NE; ++i)
				rad.jnu[id][i] += rad.jnu_line[id][i];
		return;
	}

	// output 
	std::vector<LineEscapeDiag> diag(NL);
	std::vector<std::vector<double>> thin_depth(NL, std::vector<double>(ND, 0.0));
	std::vector<std::vector<double>> esc_depth(NL, std::vector<double>(ND, 0.0));
	std::vector<std::vector<double>> heat_depth(NL, std::vector<double>(ND, 0.0));
	std::vector<std::vector<double>> beta_depth_w(NL, std::vector<double>(ND, 0.0));
	std::vector<std::vector<double>> pelec_depth_w(NL, std::vector<double>(ND, 0.0));
	std::vector<std::vector<double>> pdest_depth_w(NL, std::vector<double>(ND, 0.0));
	std::vector<std::vector<double>> y_depth_w(NL, std::vector<double>(ND, 0.0));

	// Column line opacity: absorbing cells matter even when they emit nothing.
	std::vector<std::vector<double>> kL(NL, std::vector<double>(ND, 0.0));
	for (int id = 0; id < ND; ++id)
		for (const LineRec& r : store[id])
		{
			const int L = idx[r.ip];
			kL[L][id] = r.kappaL;
			if (diag[L].ip < 0)
			{
				diag[L].ip = r.ip;
				diag[L].E_eV = r.E_eV;
				diag[L].wave = r.wave;
				diag[L].label = r.label;
				diag[L].comment = r.comment;
			}
		}

	// Cumulative optical depth for each transition:
	// tau(id -> boundary) = integral kappa_L ds through the whole column.
	std::vector<std::vector<double>> pref(NL, std::vector<double>(ND, 0.0));
	std::vector<double> tot(NL, 0.0);
	for (int L = 0; L < NL; ++L)
	{
		double c = 0.0;
		for (int id = 0; id < ND; ++id)
		{
			pref[L][id] = c;
			c += kL[L][id] * g.dr[id];
		}
		tot[L] = c;
		diag[L].tau_col_max = c;
	}

	// For each emitted line photon packet, compute the surviving fraction P.
	for (int id = 0; id < ND; ++id)
	{
		for (const LineRec& r : store[id])
		{
			const int L = idx[r.ip];
			// Emission point = cell midpoint.
			double tau_in  = pref[L][id] + 0.5 * kL[L][id] * g.dr[id];
			double tau_out = tot[L] - tau_in;
			if (tau_in  < 0.0) tau_in  = 0.0;
			if (tau_out < 0.0) tau_out = 0.0;

			// beta: two-sided line escape probability.
			double beta = 0.5 * (esc_oneside(r.redis, tau_in,  r.damp)
			                   + esc_oneside(r.redis, tau_out, r.damp));

			// Pdest: continuum absorption after the photon is trapped in the line.
			// Direct continuum attenuation is handled later by the RT solver.
			double Pdest = 0.0;
			if (beta < 1.0 && r.kappaL > 0.0)
			{
				// Continuum true absorption at h*nu_ul.
				const double conopc = rad.kabs[id][r.bin];
				if (conopc > 0.0)
				{
					const double sqrt_pi = 1.7724538509055160273;
					// eps = continuum absorption / (line-center absorption + continuum).
					const double eps = conopc / (sqrt_pi * r.kappaL + conopc);

					// Convert per-encounter eps into a trapping destruction probability.
					double dfit;
					if (r.redis == ipLY_A)
					{
						// H Lya: tau-dependent Hummer-Kunasz fit; remove Thomson
						// scattering because Pelec handles it as an escape channel.
						auto lya_side = [&](double tau) -> double
						{
							double taulog = (tau > 0.0) ?
								log10(std::min(1.0e8, sqrt_pi * tau)) : 0.0;
							double den = 0.30972 -
								std::min(0.28972, 0.03541667 * taulog);
							return eps / den;
						};
						dfit = 0.5 * (lya_side(tau_in) + lya_side(tau_out));
						const double ktot = conopc + rad.ksct[id][r.bin];
						if (ktot > 0.0)
							dfit *= conopc / ktot;
					}
					else
					{
						dfit = std::min(1.0e-3, 8.5 * eps);
						dfit /= 1.0 + dfit;
					}

					// Only trapped photons can be destroyed by continuum absorption.
					Pdest = (1.0 - beta) * dfit;
					if (Pdest < 0.0) Pdest = 0.0;
					if (Pdest > 1.0 - beta) Pdest = 1.0 - beta;
				}
			}

			// Pelec: electron scattering shifts a trapped photon out of the line core.
			double Pelec = 0.0;
			if (beta < 1.0 && r.kappaLo > 0.0)
			{
				const double kelec = rad.n_e[id] * SIGMA_THOMSON;
				Pelec = kelec / (kelec + r.kappaLo)
				      * std::max(0.0, 1.0 - beta - Pdest);
			}

			// Final surviving fraction
			// a. normal escape probability (CRD,PRD,etc.);
			// b. line scattering escape probability;
			// c. collisional de-excitation
			// d. destruction in random walk.
			double denom = beta + Pelec + r.y + Pdest;
			double P = (denom > 0.0) ?
				(beta + Pelec) * (1.0 + r.y) / denom : 1.0;
			if (P < 0.0) P = 0.0;
			if (P > 1.0) P = 1.0;

			// Assume the line photons, which destroyed by continuum absorption, 
			// are all contribute to the thermal source.
			// TODO: seperate it to thermal energy and ionization threshold
			if (denom > 0.0 && Pdest > 0.0)
				rad.line_heat[id] += r.emiss * (1.0 + r.y) * Pdest / denom;

			// These are stored for later output
			if (r.emiss > 0.0)
			{
				const double Hdest = (denom > 0.0 && Pdest > 0.0) ?
					r.emiss * (1.0 + r.y) * Pdest / denom : 0.0;

				diag[L].thin += r.emiss;
				diag[L].escaped += r.emiss * P;
				diag[L].heat_dest += Hdest;
				diag[L].beta_w += r.emiss * beta;
				diag[L].pelec_w += r.emiss * Pelec;
				diag[L].pdest_w += r.emiss * Pdest;
				diag[L].y_w += r.emiss * r.y;
				diag[L].tau_mid_w += r.emiss * g.tau_mid[id];

				thin_depth[L][id] += r.emiss;
				esc_depth[L][id] += r.emiss * P;
				heat_depth[L][id] += Hdest;
				beta_depth_w[L][id] += r.emiss * beta;
				pelec_depth_w[L][id] += r.emiss * Pelec;
				pdest_depth_w[L][id] += r.emiss * Pdest;
				y_depth_w[L][id] += r.emiss * r.y;
			}

			// Convert escaped line power to per-eV, per-sr emissivity.
			rad.jnu_line[id][r.bin] += r.emiss / (r.dE_eV * phys::four_pi) * P;
		}
	}

	// this is a diagnostic function for escape probability, used for figure 1 in paper

	// write_line_escape_diagnostics(rad, g, par, iter, diag, thin_depth, esc_depth,
	//                               heat_depth, beta_depth_w, pelec_depth_w,
	//                               pdest_depth_w, y_depth_w);

	// 5. Fold lines (B) + fluorescence (C) into jnu.
	for (int id = 0; id < ND; ++id)
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
	double E_loryd = g.E_IN_LO/phys::eV_per_Ryd;
	double E_hiryd = g.E_IN_HI/phys::eV_per_Ryd;
	snprintf(buf, sizeof(buf), "intensity %.8f range %.8f to %.8f",(rad.log_xi[id]-log10(phys::four_pi))+par.nh,E_loryd,E_hiryd);
	cdRead(buf);

	// Feed the previous RT pass's destroyed line power back into Cloudy's
	// local thermal balance. HEXTRA expects log10(erg cm^-3 s^-1).
	if (rad.line_heat[id] > 0.0 && std::isfinite(rad.line_heat[id]))
	{
		snprintf(buf, sizeof(buf), "hextra %.8f", log10(rad.line_heat[id]));
		cdRead(buf);
	}
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
	snprintf(buf, sizeof(buf), "intensity %.8f range %.8f to %.8f",(rad.log_xi[id]-log10(phys::four_pi))+par.nh,g.E_IN_LO/phys::eV_per_Ryd,g.E_IN_HI/phys::eV_per_Ryd);
	cdRead(buf);

	char fname[256];
	snprintf(fname, sizeof(fname), "\"%s%i.iron\"", par.run_hash,id);
	snprintf(buf, sizeof(buf), "save element iron %s", fname);
	cdRead(buf);
}
