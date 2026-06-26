#include "params.h"
#include "rt_grids.h"
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <sys/stat.h>

ModelParams read_params(int argc, char *argv[])
{
	// ---- default values (fill in yours) ----
	ModelParams p;
	p.nh    = 15;  // log hydrogen density [cm^-3]
	p.zeta  = 3.0;  // ionization parameter: xi = 10^zeta
	p.Gamma = -1;       // unset — required by corona model
	p.E_cut = -1;       // unset
	p.E_lo_cut = 0.1; // 0.1 keV — fixed low-energy cutoff
	p.kT_e  = -1;    // unset
	p.kT_bb = -1;    // unset
	p.taup  = -1;    // unset
	p.kT_disk = 0.35;
	strncpy(p.corona, "", sizeof(p.corona));
	p.frac      = -1;   // flux ratio: F_corona/F_disk. If frac <= 0, no I_disk (corona only, default)
	p.incidence = 0.7071067811865476;  // cos(45 deg)
	p.test_rt   = false;
	strncpy(p.test_mode, "none", sizeof(p.test_mode));
	p.tau_slab  = 0.5;   // total vertical Thomson optical depth (compps mode)
	p.angsca = true;
	p.E_rt_lo   = 10.0;      // 0.01 keV in eV
	p.E_rt_hi   = 1000.0e3;   // 1000 keV in eV
	p.maxiter   = 500;
	p.Afe       = 1.0;     // Fe/Fe_solar (linear), 1.0 = solar
	p.run_hash[0] = '\0';
	p.run_dir[0]  = '\0';
	// Approximated kernel has some issue; we need fix. (ktype=0)
	// It never match the normalization condition.
	p.ktype = 1;	

	// ---- command-line overrides ----
	for (int i = 1; i < argc; ++i)
	{
		if (strcmp(argv[i], "-nh") == 0 && i+1 < argc)
			p.nh = atof(argv[++i]);
		else if (strcmp(argv[i], "-zeta") == 0 && i+1 < argc)
			p.zeta = atof(argv[++i]);
		else if (strcmp(argv[i], "-Gamma") == 0 && i+1 < argc)
			p.Gamma = atof(argv[++i]);
		else if (strcmp(argv[i], "-Ecut") == 0 && i+1 < argc)
			p.E_cut = atof(argv[++i]);
		else if (strcmp(argv[i], "-E_low_cut") == 0 && i+1 < argc)
			p.E_lo_cut = atof(argv[++i]);
		else if (strcmp(argv[i], "-frac") == 0 && i+1 < argc)
			p.frac = atof(argv[++i]);
		else if (strcmp(argv[i], "-incidence") == 0 && i+1 < argc)
			p.incidence = atof(argv[++i]);
		else if (strcmp(argv[i], "-test_rt") == 0)
		{
			p.test_rt = true;
			// required mode token: "-test_rt compps" (or test_avg)
			if (i+1 < argc && argv[i+1][0] != '-')
				strncpy(p.test_mode, argv[++i], sizeof(p.test_mode) - 1);
		}
		else if (strcmp(argv[i], "-tau") == 0 && i+1 < argc)
			p.tau_slab = atof(argv[++i]);
		else if (strcmp(argv[i], "-kT_disk") == 0 && i+1 < argc)
			p.kT_disk = atof(argv[++i]);
		else if (strcmp(argv[i], "-corona") == 0 && i+1 < argc)
			strncpy(p.corona, argv[++i], sizeof(p.corona) - 1);
		else if (strcmp(argv[i], "-kT_e") == 0 && i+1 < argc)
			p.kT_e = atof(argv[++i]);
		else if (strcmp(argv[i], "-kT_bb") == 0 && i+1 < argc)
			p.kT_bb = atof(argv[++i]);
		else if (strcmp(argv[i], "-taup") == 0 && i+1 < argc)
			p.taup = atof(argv[++i]);
		else if (strcmp(argv[i], "-Afe") == 0 && i+1 < argc)
			p.Afe = atof(argv[++i]);
		else if (strcmp(argv[i], "-angsca") == 0 && i+1 < argc)
		{
			// Angular-scattering kernel selector:
			//   1/true/yes → angle-dependent KernelCache
			//   0/false/no → angle-mean   avgKernelCache
			const char* v = argv[++i];
			p.angsca = (strcmp(v, "1") == 0 || strcmp(v, "true") == 0
			         || strcmp(v, "yes") == 0);
		}
	}

	// --- Validate corona model and its required parameters ---
	if (p.corona[0] == '\0')
	{
		fprintf(stderr, "Error: -corona is required. Options: powerlaw, cutoffpl, nthcomp, comptt, blackbody\n");
		fprintf(stderr, "  -corona powerlaw   requires: -Gamma\n");
		fprintf(stderr, "  -corona cutoffpl   requires: -Gamma -Ecut\n");
		fprintf(stderr, "  -corona nthcomp    requires: -Gamma -kT_e -kT_bb\n");
		fprintf(stderr, "  -corona comptt     requires: -kT_e -kT_bb -taup\n");
		fprintf(stderr, "  -corona blackbody  requires: -kT_bb\n");
		exit(1);
	}
	else if (strcmp(p.corona, "powerlaw") == 0)
	{
		if (p.Gamma < 0) {
			fprintf(stderr, "Error: corona=powerlaw requires -Gamma\n"); exit(1);
		}
		printf("Corona:     powerlaw  Gamma=%.4f\n", p.Gamma);
	}
	else if (strcmp(p.corona, "cutoffpl") == 0)
	{
		if (p.Gamma < 0 || p.E_cut < 0) {
			fprintf(stderr, "Error: corona=cutoffpl requires -Gamma -Ecut\n"); exit(1);
		}
		printf("Corona:     cutoffpl  Gamma=%.4f  E_cut=%.2f keV\n", p.Gamma, p.E_cut);
	}
	else if (strcmp(p.corona, "nthcomp") == 0)
	{
		if (p.Gamma < 0 || p.kT_e < 0 || p.kT_bb < 0) {
			fprintf(stderr, "Error: corona=nthcomp requires -Gamma -kT_e -kT_bb\n"); exit(1);
		}
		printf("Corona:     nthcomp  Gamma=%.4f  kT_e=%.2f keV  kT_bb=%.4f keV\n",
		       p.Gamma, p.kT_e, p.kT_bb);
	}
	else if (strcmp(p.corona, "comptt") == 0)
	{
		if (p.kT_e < 0 || p.kT_bb < 0 || p.taup < 0) {
			fprintf(stderr, "Error: corona=comptt requires -kT_e -kT_bb -taup\n"); exit(1);
		}
		printf("Corona:     comptt  kT_e=%.2f keV  kT_bb=%.4f keV  taup=%.2f\n",
		       p.kT_e, p.kT_bb, p.taup);
	}
	else if (strcmp(p.corona, "blackbody") == 0)
	{
		if (p.kT_bb < 0) {
			fprintf(stderr, "Error: corona=blackbody requires -kT_bb\n"); exit(1);
		}
		printf("Corona:     blackbody  kT_bb=%.4f keV\n", p.kT_bb);
	}
	else
	{
		fprintf(stderr, "Error: unknown corona model '%s'\n", p.corona);
		fprintf(stderr, "  Available: powerlaw, cutoffpl, nthcomp, comptt, blackbody\n");
		exit(1);
	}

	printf("Kernel:     %s\n",
	       p.angsca ? "angle-dependent (KernelCache)"
	                : "angle-mean (avgKernelCache)");

	double xi = pow(10.0, p.zeta);
	double nH = pow(10.0, p.nh);
	double Fx = xi * nH / 4.0 / M_PI ;
	printf("Parameters: nh=%.4f  zeta=%.4f (xi=%.4e)  frac=%.4f  incidence=%.4f  Afe=%.2f\n",
	       p.nh, p.zeta, xi, p.frac, p.incidence, p.Afe);
	printf("Derived:    nH=%.4e  Fx=%.4e erg/cm^2/s\n", nH, Fx);

	// --- Compute deterministic run hash from physics params ---
	{
		char buf[512];
		snprintf(buf, sizeof(buf),
			"%s|mode=%s|nh=%.6g|zeta=%.6g|frac=%.6g|inc=%.6g|Afe=%.6g|"
			"Gamma=%.6g|Ecut=%.6g|Elo=%.6g|kTe=%.6g|kTbb=%.6g|"
			"taup=%.6g|kTd=%.6g|test=%d|tau=%.6g|ang=%d",
			p.corona, p.test_mode, p.nh, p.zeta, p.frac, p.incidence, p.Afe,
			p.Gamma, p.E_cut, p.E_lo_cut, p.kT_e, p.kT_bb,
			p.taup, p.kT_disk, (int)p.test_rt, p.tau_slab,
			(int)p.angsca);
		// FNV-1a 32-bit hash
		unsigned int h = 2166136261u;
		for (const char* c = buf; *c; ++c)
		{
			h ^= (unsigned char)*c;
			h *= 16777619u;
		}
		snprintf(p.run_hash, sizeof(p.run_hash), "%08x", h);
		snprintf(p.run_dir, sizeof(p.run_dir), "results/%s", p.run_hash);
		printf("Run hash:   %s  →  %s/\n", p.run_hash, p.run_dir);
	}

	// --- Create run directory and write params.json ---
	{
		mkdir("results", 0755);
		mkdir(p.run_dir, 0755);

		char pjson[512];
		snprintf(pjson, sizeof(pjson), "%s/params.json", p.run_dir);
		FILE* fp = fopen(pjson, "w");
		if (fp)
		{
			char timestamp[32];
			time_t now = time(nullptr);
			strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S",
			         localtime(&now));

			fprintf(fp, "{\n");
			fprintf(fp, "  \"hash\": \"%s\",\n", p.run_hash);
			fprintf(fp, "  \"time\": \"%s\",\n", timestamp);
			fprintf(fp, "  \"corona\": \"%s\",\n", p.corona);
			fprintf(fp, "  \"nh\": %.6g,\n", p.nh);
			fprintf(fp, "  \"zeta\": %.6g,\n", p.zeta);
			fprintf(fp, "  \"frac\": %.6g,\n", p.frac);
			fprintf(fp, "  \"incidence\": %.6g,\n", p.incidence);
			fprintf(fp, "  \"Afe\": %.6g,\n", p.Afe);
			fprintf(fp, "  \"Gamma\": %.6g,\n", p.Gamma);
			fprintf(fp, "  \"E_cut\": %.6g,\n", p.E_cut);
			fprintf(fp, "  \"E_lo_cut\": %.6g,\n", p.E_lo_cut);
			fprintf(fp, "  \"kT_e\": %.6g,\n", p.kT_e);
			fprintf(fp, "  \"kT_bb\": %.6g,\n", p.kT_bb);
			fprintf(fp, "  \"taup\": %.6g,\n", p.taup);
			fprintf(fp, "  \"kT_disk\": %.6g,\n", p.kT_disk);
			fprintf(fp, "  \"test_rt\": %s,\n", p.test_rt ? "true" : "false");
			fprintf(fp, "  \"test_mode\": \"%s\",\n", p.test_mode);
			fprintf(fp, "  \"tau_slab\": %.6g,\n", p.tau_slab);
			fprintf(fp, "  \"maxiter\": %d,\n", p.maxiter);
			fprintf(fp, "  \"E_rt_lo\": %.6g,\n", p.E_rt_lo);
			fprintf(fp, "  \"E_rt_hi\": %.6g,\n", p.E_rt_hi);
			fprintf(fp, "  \"angsca\": %s,\n", p.angsca ? "true" : "false");
			fprintf(fp, "  \"ktype\": %d\n", p.ktype);
			fprintf(fp, "}\n");
			fclose(fp);
		}
	}

	return p;
}

void snap_incidence(ModelParams& par, const RTGrids& g)
{
	// snap incidence to nearest GL node (negative = downward)
	double mu_target = -fabs(par.incidence);
	int i_mu = 0;
	double best = fabs(mu_target - g.mu[0]);
	for (int j = 1; j < RTGrids::NA; ++j)
	{
		double diff = fabs(mu_target - g.mu[j]);
		if (diff < best) { best = diff; i_mu = j; }
	}
	par.incidence   = g.mu[i_mu];
	par.i_incidence = i_mu;

	printf("Incidence:  cos(theta)=%.6f  (g.mu[%d])\n", par.incidence, par.i_incidence);
}
