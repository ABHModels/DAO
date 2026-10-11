#include "test_rt.h"
#include "rt_grids.h"
#include "constants.h"
#include "compton_rt.h"
#include "corona_models.h"   // blackbody()

#include <cmath>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include "run_log.h"
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>

// ============================================================
// Fill data for RT test mode: isothermal pure-scattering slab.
// ============================================================
static void fill_test_data(RadField& rad, const RTGrids& g,
                           double nh, double T_slab)
{
	const double nH  = pow(10.0, nh);
	const double n_e = phys::reference_electrons_per_hydrogen * nH;

	// The shared scattering operator supplies opacity when the slab is solved.

	for (int id = 0; id < g.ND_MID; ++id)
	{
		rad.T_K[id] = T_slab;
		rad.n_e[id] = n_e;

		for (int ie = 0; ie < g.NE; ++ie)
		{
			rad.kabs[id][ie] = 0.0;
			rad.jnu[id][ie]  = 0.0;
		}


		rad.heating[id] = 0.0;
		rad.line_heat[id] = 0.0;
		rad.cooling[id] = 0.0;
	}

	double kb_eV = phys::k_B / phys::eV_to_erg;
	fprintf(stdout, "Test mode: pure scattering atmosphere\n");
	fprintf(stdout, "  T=%.1e K  kT=%.2f eV  n_e=%.3e cm^-3\n",
		T_slab, kb_eV * T_slab, n_e);
}

// Prescribed-temperature tests share production's RT and scattering operator.
// Hot electrons exchange energy with photons, so Fout need not equal Fin.
// Photon number and the signed boundary/volume energy identity must still hold.
template<class Cache>
static void solve_test_slab(RadField& rad,const RTGrids& g,
                            const ModelParams& par,const Cache& cache)
{
    const auto ops=make_scattering_column(rad,g,cache);
    compton_rt_solve(rad,g,par,ops);
    const auto w=dao_thermal::energy_weights(g.NE,g.ene);
    double fin=0,fout=0,nin=0,nout=0,volume=0;
    for(int m=0;m<g.NA;++m) for(int e=0;e<g.NE;++e) {
        const double weight=0.5*phys::four_pi*g.wt[m]*std::abs(g.mu[m])*w[e];
        const double in=g.mu[m]<0 ? rad.Inu_top[m][e] : rad.Inu_bottom[m][e];
        const double out=g.mu[m]<0 ? rad.Inu_bottom[m][e] : rad.Inu_top[m][e];
        fin+=weight*in;fout+=weight*out;
        nin+=weight*in/(g.ene[e]*phys::eV_to_erg);
        nout+=weight*out/(g.ene[e]*phys::eV_to_erg);
    }
    std::vector<double> intensity(g.NA*g.NE);
    for(int d=0;d<g.ND_MID;++d) {
        for(int m=0;m<g.NA;++m)
            std::copy(rad.Inu[d][m],rad.Inu[d][m]+g.NE,intensity.data()+m*g.NE);
        volume+=dao_thermal::budget(g.NE,w.data(),rad.J0[d],rad.kabs[d],rad.jnu[d],ops[d],
                                    par.angsca ? intensity.data() : nullptr).residual()*g.dr[d];
    }
    const double photons=nout/nin-1,identity=(fout-fin+volume)/fin;
    if(!std::isfinite(photons) || !std::isfinite(identity) ||
       std::abs(photons)>1e-6 || std::abs(identity)>1e-6)
        throw std::runtime_error("RT benchmark: photon or boundary/volume energy conservation failed");
    std::ofstream file(std::string(par.run_dir)+"/test_rt_budget.json");
    if(!file) throw std::runtime_error("Cannot write RT benchmark budget");
    file<<std::setprecision(16)<<"{\n  \"incoming_flux\": "<<fin
        <<",\n  \"outgoing_flux\": "<<fout<<",\n  \"volume_residual\": "<<volume
        <<",\n  \"incoming_photons\": "<<nin<<",\n  \"outgoing_photons\": "<<nout
        <<",\n  \"relative_photon_error\": "<<photons
        <<",\n  \"energy_identity_error\": "<<identity<<"\n}\n";
    dao_log::info("RT check:   photon error %.3e | energy identity %.3e\n",photons,identity);
}

// ============================================================
// compps benchmark mode
//
// Isothermal, pure-scattering plane-parallel slab illuminated
// from the BOTTOM by an isotropic blackbody seed — the exact
// configuration of Xspec's compps with:
//   geom=1 (slab), cov_frac=1 (clean slab), rel_refl=0 (no reflection),
//   Maxwellian electrons at kTe, seed blackbody at kTbb.
//
// Parameter mapping (compps  ->  our flags):
//   kTe   [keV]  ->  -kT_e   (uniform electron temperature T_K)
//   kTbb  [keV]  ->  -kT_bb  (bottom blackbody seed)
//   tau          ->  -tau    (total Thomson depth, set on the grid in main)
//   cosIncl      ->  -incidence (emergent intensity reported at GL nodes)
// ============================================================
static void save_emergent_compps(const RadField& rad, const RTGrids& g,
                                  const ModelParams& par, double Te_keV)
{
	// Save into the run directory, same logic as save_results() in the
	// production path: results/<hash>/.  Column layout matches the
	// production emergent_iter file so existing readers work.
	const char* dir = par.run_dir[0] ? par.run_dir : "results";
	mkdir("results", 0755);
	if (par.run_dir[0])
		mkdir(dir, 0755);

	char fname[256];
	snprintf(fname, sizeof(fname), "%s/emergent_compps.dat", dir);

	int fd = open(fname, O_WRONLY | O_CREAT | O_TRUNC, 0644);
	FILE* fp = fdopen(fd, "w");

	fprintf(fp, "# compps benchmark slab (isothermal pure scattering)\n");
	fprintf(fp, "# kTe=%.3f keV  kTbb=%.3f keV  tau_T=%.3f  geom=slab\n",
	        Te_keV, par.kT_bb, par.tau_slab);
	fprintf(fp, "# Col 1: E [eV]\n");
	fprintf(fp, "# Col 2: F_corona (incident normal flux) [erg cm^-2 s^-1 eV^-1]\n");
	fprintf(fp, "# Col 3: I_disk   (bottom blackbody seed, incident)\n");
	for (int nm = 0; nm < g.NA; ++nm)
		fprintf(fp, "# Col %d: I_emergent(mu=%.6f)\n", nm + 4, g.mu[nm]);

	for (int ie = 0; ie < g.NE; ++ie)
	{
		fprintf(fp, "%.6e  %.6e  %.6e",
		        g.ene[ie], rad.illum.I_corona[ie], rad.illum.I_disk[ie]);
		for (int nm = 0; nm < g.NA; ++nm)
			fprintf(fp, "  %.6e", rad.Inu_top[nm][ie]);
		fprintf(fp, "\n");
	}

	fclose(fp);
	fprintf(stdout, "Saved: %s\n", fname);
}

static void save_emergent_tavg(const RadField& rad, const RTGrids& g,
                                  const ModelParams& par, double Te_keV)
{
	// Save into the run directory, same logic as save_results() in the
	// production path: results/<hash>/.  Column layout matches the
	// production emergent_iter file so existing readers work.
	const char* dir = par.run_dir[0] ? par.run_dir : "results";
	mkdir("results", 0755);
	if (par.run_dir[0])
		mkdir(dir, 0755);

	char fname[256];
	snprintf(fname, sizeof(fname), "%s/emergent_tavg.dat", dir);

	int fd = open(fname, O_WRONLY | O_CREAT | O_TRUNC, 0644);
	FILE* fp = fdopen(fd, "w");

	fprintf(fp, "# compps benchmark slab (isothermal pure scattering)\n");
	fprintf(fp, "# kTe=%.3f keV  kTbb=%.3f keV  tau_T=%.3f  geom=slab\n",
	        Te_keV, par.kT_bb, par.tau_slab);
	fprintf(fp, "# Col 1: E [eV]\n");
	fprintf(fp, "# Col 2: F_corona (incident normal flux) [erg cm^-2 s^-1 eV^-1]\n");
	fprintf(fp, "# Col 3: I_disk   (bottom blackbody seed, incident)\n");
	for (int nm = 0; nm < g.NA; ++nm)
		fprintf(fp, "# Col %d: I_emergent(mu=%.6f)\n", nm + 4, g.mu[nm]);

	for (int ie = 0; ie < g.NE; ++ie)
	{
		fprintf(fp, "%.6e  %.6e  %.6e",
		        g.ene[ie], rad.illum.I_corona[ie], rad.illum.I_disk[ie]);
		for (int nm = 0; nm < g.NA; ++nm)
			fprintf(fp, "  %.6e", rad.Inu_top[nm][ie]);
		fprintf(fp, "\n");
	}

	fclose(fp);
	fprintf(stdout, "Saved: %s\n", fname);
}

template<class Cache>
static void run_test_rt_compps(RadField& rad, const RTGrids& g,
                               ModelParams& par, Cache& kcache)
{
	// kTe must be supplied (the slab electron temperature) and, for a
	// Maxwellian, compps itself requires kTe > 1 keV.
	if (par.kT_e <= 0.0)
	{
		fprintf(stderr,
		        "Error: -test_rt compps requires -kT_e (slab Te in keV)\n");
		exit(1);
	}

	// kTe [keV] -> uniform slab temperature [K]
	const double kB_eV = phys::k_B / phys::eV_to_erg;   // [eV/K]
	const double Te_K  = par.kT_e * 1.0e3 / kB_eV;

	fprintf(stdout, "\n=== RT TEST MODE: compps benchmark ===\n");
	fprintf(stdout, "  kTe=%.3f keV (T=%.3e K)  kTbb=%.3f keV  tau_T=%.3f\n",
	        par.kT_e, Te_K, par.kT_bb, par.tau_slab);

	// Isothermal pure-scattering slab at Te.
	fill_test_data(rad, g, par.nh, Te_K);

	// Illumination: no top source; isotropic blackbody seed at kTbb
	// enters from the bottom with specific intensity ill_bot = 2*I_disk.
	blackbody(rad.illum.I_disk, g, par.kT_bb * 1.0e3);
	for (int ie = 0; ie < g.NE; ++ie)
		rad.illum.I_corona[ie] = 0.0;

	// Prescribed bottom-seed mean intensity for the compPS shape benchmark.
	// Retain its reference amplitude; this is not coronal illumination.
	double raw_disk = 0.0;
	for (int ie = 0; ie < g.NE - 1; ++ie)
	{
		double dE = g.ene[ie + 1] - g.ene[ie];
		if (g.ene[ie] > g.E_IN_LO && g.ene[ie] < g.E_IN_HI)
			raw_disk += 0.5 * (rad.illum.I_disk[ie] + rad.illum.I_disk[ie + 1]) * dE;
	}
	double xi = pow(10.0, par.zeta);
	double nH = pow(10.0, par.nh);
	double J_seed = xi * nH / pow(phys::four_pi, 2);
	double scale_disk = J_seed / raw_disk;
	for (int ie = 0; ie < g.NE; ++ie)
		rad.illum.I_disk[ie] *= scale_disk;

	fprintf(stdout, "  Seed mean intensity: J_seed=%.4e  raw_disk=%.4e\n", J_seed, raw_disk);

	solve_test_slab(rad, g, par, kcache);
	rad.compute_moments();
	save_emergent_compps(rad, g, par, par.kT_e);
}

// ============================================================
// test_avg mode
//
// Isothermal, pure-scattering slab illuminated from the top by the production
// corona spectrum. Uses the same fixed incident-flux normalization and boundary
// conversion as production; no bottom source.
// ============================================================
template<class Cache>
static void run_test_rt_avg(RadField& rad, const RTGrids& g,
                            ModelParams& par, Cache& kcache)
{
	if (par.kT_e <= 0.0)
	{
		fprintf(stderr,
		        "Error: -test_rt test_avg requires -kT_e (slab Te in keV)\n");
		exit(1);
	}

	const double kB_eV = phys::k_B / phys::eV_to_erg;   // [eV/K]
	const double Te_K  = par.kT_e * 1.0e3 / kB_eV;

	fprintf(stdout, "\n=== RT TEST MODE: test_avg (top pencil beam) ===\n");
	fprintf(stdout, "  corona=%s  kTe=%.3f keV (T=%.3e K)  tau_T=%.3f\n",
	        par.corona, par.kT_e, Te_K, par.tau_slab);

	// Isothermal pure-scattering slab at Te (same as compps).
	fill_test_data(rad, g, par.nh, Te_K);

	// Use the production illumination path, with the bottom source disabled.
	ModelParams illumination = par;
	illumination.frac = -1;
	rad.illum.compute(illumination);

	solve_test_slab(rad, g, par, kcache);
	rad.compute_moments();
	save_emergent_tavg(rad, g, par, par.kT_e);
}

// ============================================================
// run_test_rt
// ============================================================
template<class Cache>
void run_test_rt(RadField& rad, const RTGrids& g, ModelParams& par,
                 Cache& kcache)
{	
	if (strcmp(par.test_mode, "compps") == 0)
	{
		run_test_rt_compps(rad, g, par, kcache);
		return;
	}
	if (strcmp(par.test_mode, "test_avg") == 0)
	{
		run_test_rt_avg(rad, g, par, kcache);
		return;
	}
	fprintf(stderr,
	        "Error: unknown test mode '%s' (valid: compps, test_avg)\n",
	        par.test_mode);
	exit(1);
}

// Explicit instantiations for both kernel-cache types.
template void run_test_rt<KernelCache>(
	RadField&, const RTGrids&, ModelParams&, KernelCache&);
template void run_test_rt<avgKernelCache>(
	RadField&, const RTGrids&, ModelParams&, avgKernelCache&);
