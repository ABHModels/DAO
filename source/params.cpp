#include "params.h"
#include "rt_grids.h"
#include "run_log.h"
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <sys/stat.h>

// Build a one-line human-readable description of the run from its physics
// params (mirrors the per-model console formatting in read_params).
static void build_run_summary(const ModelParams& p, char* out, size_t n)
{
	char model[160];
	if (p.test_rt)
		snprintf(model, sizeof(model), "%s test | tau=%.4g | %s",
		         p.test_mode, p.tau_slab, p.corona);
	else if (strcmp(p.corona, "powerlaw") == 0)
		snprintf(model, sizeof(model), "powerlaw | Gamma=%.4g", p.Gamma);
	else if (strcmp(p.corona, "cutoffpl") == 0)
		snprintf(model, sizeof(model), "cutoffpl | Gamma=%.4g Ecut=%.4g keV",
		         p.Gamma, p.E_cut);
	else if (strcmp(p.corona, "nthcomp") == 0)
		snprintf(model, sizeof(model),
		         "nthcomp | Gamma=%.4g kTe=%.4g kTbb=%.4g keV",
		         p.Gamma, p.kT_e, p.kT_bb);
	else if (strcmp(p.corona, "comptt") == 0)
		snprintf(model, sizeof(model),
		         "comptt | kTe=%.4g kTbb=%.4g keV taup=%.4g",
		         p.kT_e, p.kT_bb, p.taup);
	else if (strcmp(p.corona, "blackbody") == 0)
		snprintf(model, sizeof(model), "blackbody | kTbb=%.4g keV", p.kT_bb);
	else
		snprintf(model, sizeof(model), "%s", p.corona);

	snprintf(out, n, "%s | nh=%.4g zeta=%.4g frac=%.4g | %s | reflionx flux normalization | depth=%s",
	         model, p.nh, p.zeta, p.frac,
	         p.angsca ? "angle-dependent" : "angle-mean",
	         p.log_depth ? "log" : "tanh");
}

// Escape all JSON control characters as well as quotes and backslashes.
static std::string json_escape(const std::string& text)
{
    std::string out;
    for (unsigned char c : text) {
        if (c == '"' || c == '\\') { out += '\\'; out += char(c); }
        else if (c < 0x20) {
            char escaped[7];
            snprintf(escaped, sizeof(escaped), "\\u%04x", unsigned(c));
            out += escaped;
        } else out += char(c);
    }
    return out;
}

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
	p.frac      = -1;   // corona/disk amplitude ratio; <=0: corona only
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
		else if (strcmp(argv[i], "-log_depth") == 0)
		{
			if (i+1 >= argc) {
				fprintf(stderr, "Error: -log_depth requires 0/1, false/true, or no/yes.\n");
				exit(1);
			}
			const char* value = argv[++i];
			if (strcmp(value, "1") == 0 || strcmp(value, "true") == 0 || strcmp(value, "yes") == 0)
				p.log_depth = true;
			else if (strcmp(value, "0") == 0 || strcmp(value, "false") == 0 || strcmp(value, "no") == 0)
				p.log_depth = false;
			else {
				fprintf(stderr, "Error: -log_depth requires 0/1, false/true, or no/yes.\n");
				exit(1);
			}
		}
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
		else if (strcmp(argv[i], "-O") == 0 || strcmp(argv[i], "-Fe") == 0)
		{
			const char* flag = argv[i];
			bool enabled = true;
			if (i+1 < argc && argv[i+1][0] != '-') {
				const char* value = argv[++i];
				if (strcmp(value,"1")==0 || strcmp(value,"true")==0 || strcmp(value,"yes")==0)
					enabled = true;
				else if (strcmp(value,"0")==0 || strcmp(value,"false")==0 || strcmp(value,"no")==0)
					enabled = false;
				else {
					fprintf(stderr,"Error: %s accepts an optional 0/1, false/true, or no/yes.\n",flag);
					exit(1);
				}
			}
			if (strcmp(flag,"-O")==0) p.save_oxygen = enabled;
			else p.save_iron = enabled;
		}
		else if (strcmp(argv[i], "-verbose") == 0)
		{
			p.verbose = true;
			if (i+1 < argc && argv[i+1][0] != '-') {
				const char* value = argv[++i];
				if (strcmp(value,"1")==0 || strcmp(value,"true")==0 || strcmp(value,"yes")==0)
					p.verbose = true;
				else if (strcmp(value,"0")==0 || strcmp(value,"false")==0 || strcmp(value,"no")==0)
					p.verbose = false;
				else {
					fprintf(stderr,"Error: -verbose accepts an optional 0/1, false/true, or no/yes.\n");
					exit(1);
				}
			}
		}
		else if (strcmp(argv[i], "-label") == 0)
		{
			if (i + 1 >= argc) { fprintf(stderr, "Error: -label requires text.\n"); exit(1); }
			p.label = argv[++i];
		}
		else if (strcmp(argv[i], "-angsca") == 0)
		{
			// Angular-scattering kernel selector:
			//   1/true/yes → angle-dependent KernelCache
			//   0/false/no → angle-mean   avgKernelCache
			if (i+1 >= argc) {
				fprintf(stderr, "Error: -angsca requires 0/1, false/true, or no/yes.\n"); exit(1);
			}
			const char* v = argv[++i];
			if (strcmp(v, "1") == 0 || strcmp(v, "true") == 0 || strcmp(v, "yes") == 0)
				p.angsca = true;
			else if (strcmp(v, "0") == 0 || strcmp(v, "false") == 0 || strcmp(v, "no") == 0)
				p.angsca = false;
			else {
				fprintf(stderr, "Error: -angsca requires 0/1, false/true, or no/yes.\n"); exit(1);
			}
		}
		else {
			fprintf(stderr, "Error: unknown or incomplete option '%s'.\n", argv[i]);
			exit(1);
		}
	}

    if(!p.test_rt) {
        // Production actually selects Cloudy's mesh using RTGrids::E_LO/E_HI.
        // Report that window for the new mode instead of the old test-grid label.
        p.E_rt_lo=RTGrids::E_LO;
        p.E_rt_hi=RTGrids::E_HI;
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
	}
	else if (strcmp(p.corona, "cutoffpl") == 0)
	{
		if (p.Gamma < 0 || p.E_cut < 0) {
			fprintf(stderr, "Error: corona=cutoffpl requires -Gamma -Ecut\n"); exit(1);
		}
	}
	else if (strcmp(p.corona, "nthcomp") == 0)
	{
		if (p.Gamma < 0 || p.kT_e < 0 || p.kT_bb < 0) {
			fprintf(stderr, "Error: corona=nthcomp requires -Gamma -kT_e -kT_bb\n"); exit(1);
		}
	}
	else if (strcmp(p.corona, "comptt") == 0)
	{
		if (p.kT_e < 0 || p.kT_bb < 0 || p.taup < 0) {
			fprintf(stderr, "Error: corona=comptt requires -kT_e -kT_bb -taup\n"); exit(1);
		}
	}
	else if (strcmp(p.corona, "blackbody") == 0)
	{
		if (p.kT_bb < 0) {
			fprintf(stderr, "Error: corona=blackbody requires -kT_bb\n"); exit(1);
		}
	}
	else
	{
		fprintf(stderr, "Error: unknown corona model '%s'\n", p.corona);
		fprintf(stderr, "  Available: powerlaw, cutoffpl, nthcomp, comptt, blackbody\n");
		exit(1);
	}

	double xi = pow(10.0, p.zeta);
	if (p.test_rt && (p.save_oxygen || p.save_iron)) {
		fprintf(stderr,"Error: -O/-Fe require production mode; test_rt does not run Cloudy.\n");
		exit(1);
	}
	double nH = pow(10.0, p.nh);
	double Fx = xi * nH / 4.0 / M_PI ;

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
		// Separate depth resolutions as well as normalization and grid mappings.
		std::string hash_input(buf);
		// Preserve the identity of existing reflionx-normalized production runs.
		hash_input += "|reflionx_norm=1";
		if (!p.log_depth) hash_input += "|log_depth=0";
		// The scattering prescription changes the physics and must not overwrite
		// results made with the former fixed 1.21*nH density.
		if (!p.test_rt) hash_input += "|rt_thermal_balance=v1|scattering_density=cloudy_ne|temperature_iteration=local_cell_response_v3|initialization=boundary_intensity_v1";
		if (p.test_rt) hash_input += "|transfer=constant_cell_v1|scattering=conservative_v1|angles=" + std::to_string(RTGrids::NA);
		if (p.test_rt && strcmp(p.test_mode,"test_avg")==0) hash_input += "|seed_normalization=flux_v1";
		if (!p.test_rt && RTGrids::NA!=8) hash_input += "|angles=" + std::to_string(RTGrids::NA);
		if (!p.test_rt) hash_input += "|column_energy_tolerance=" + std::to_string(ModelParams::thermal_column_tolerance);
		hash_input += "|depth_cells=" + std::to_string(RTGrids::ND_MID_DEFAULT);
		unsigned int h = 2166136261u;
		for (unsigned char c : hash_input)
		{
			h ^= c;
			h *= 16777619u;
		}
		snprintf(p.run_hash, sizeof(p.run_hash), "%08x", h);
		snprintf(p.run_dir, sizeof(p.run_dir), "results/%s", p.run_hash);
	}

	char summary[256];
	build_run_summary(p, summary, sizeof(summary));
    if(!p.test_rt) {
        std::string label=std::string(summary)+" | RT thermal balance";
        snprintf(summary,sizeof(summary),"%s",label.c_str());
    }
	char timestamp[32];
	time_t now = time(nullptr);
	strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", localtime(&now));

	// --- Create run directory and write params.json ---
	{
		mkdir("results", 0755);
		mkdir(p.run_dir, 0755);
		dao_log::open(std::string(p.run_dir)+"/run.log",p.verbose);
		dao_log::detail("\n=== Run started %s ===\n",timestamp);
		dao_log::info("Run hash:   %s  →  %s/\n",p.run_hash,p.run_dir);
		dao_log::info("Run:        %s\n",summary);
		if(!p.label.empty()) dao_log::info("Label:      %s\n",p.label.c_str());
		dao_log::info("Normalization: %s | Afe=%.4g\n","reflionx flux",p.Afe);
		dao_log::info("Log:        %s/run.log%s\n",p.run_dir,p.verbose ? " (verbose)" : " (use -verbose for details)");
		dao_log::detail("Parameters: nh=%.4f zeta=%.4f (xi=%.4e) frac=%.4f incidence=%.6f\n",
		                p.nh,p.zeta,xi,p.frac,p.incidence);
		dao_log::detail("Derived: nH=%.4e illumination normalization=%.4e erg/cm^2/s\n",nH,Fx);

		char pjson[512];
		snprintf(pjson, sizeof(pjson), "%s/params.json", p.run_dir);
		FILE* fp = fopen(pjson, "w");
		if (fp)
		{

			fprintf(fp, "{\n");
			fprintf(fp, "  \"hash\": \"%s\",\n", p.run_hash);
			fprintf(fp, "  \"time\": \"%s\",\n", timestamp);
			fprintf(fp, "  \"label\": \"%s\",\n", json_escape(p.label).c_str());
			fprintf(fp, "  \"verbose\": %s,\n", p.verbose ? "true" : "false");
			fprintf(fp, "  \"save_oxygen\": %s,\n", p.save_oxygen ? "true" : "false");
			fprintf(fp, "  \"save_iron\": %s,\n", p.save_iron ? "true" : "false");
			fprintf(fp, "  \"corona\": \"%s\",\n", p.corona);
			fprintf(fp, "  \"nh\": %.6g,\n", p.nh);
			fprintf(fp, "  \"zeta\": %.6g,\n", p.zeta);
			fprintf(fp, "  \"frac\": %.6g,\n", p.frac);
			fprintf(fp, "  \"incidence\": %.6g,\n", p.incidence);
			fprintf(fp, "  \"illumination_normalization\": \"reflionx\",\n");
			fprintf(fp, "  \"log_depth\": %s,\n", p.log_depth ? "true" : "false");
			fprintf(fp, "  \"depth_cells\": %d,\n", RTGrids::ND_MID_DEFAULT);
            fprintf(fp, "  \"transfer_solver\": \"constant_cell_v1\",\n");
            fprintf(fp, "  \"angular_points\": %d,\n", RTGrids::NA);
            fprintf(fp, "  \"intensity_location\": \"cell-volume-average\",\n");
            // Read-only metadata for old result readers, not a solver switch.
            fprintf(fp, "  \"rt_thermal_balance\": %s,\n", !p.test_rt ? "true" : "false");
            if(!p.test_rt) {
                fprintf(fp, "  \"solver\": \"rt_energy_balance_v1\",\n");
                fprintf(fp, "  \"scattering_density\": \"cloudy_free_electrons\",\n");
                fprintf(fp, "  \"temperature_iteration\": \"local_cell_response_v3\",\n");
                fprintf(fp, "  \"column_energy_tolerance\": %.9g,\n", ModelParams::thermal_column_tolerance);
                fprintf(fp, "  \"initialization\": \"boundary_intensity_v1\",\n");
            }
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
		// --- Human-readable run summary: RUN.txt ---
		char rtxt[512];
		snprintf(rtxt, sizeof(rtxt), "%s/RUN.txt", p.run_dir);
		FILE* rf = fopen(rtxt, "w");
		if (!rf) { perror(rtxt); exit(1); }
		{
			fprintf(rf, "Run:     %s\n", summary);
			fprintf(rf, "Label:   %s\n", !p.label.empty() ? p.label.c_str() : "(none)");
			fprintf(rf, "Time:    %s\n", timestamp);
			fprintf(rf, "Hash:    %s\n", p.run_hash);
			fprintf(rf, "Final ion fractions: O=%s Fe=%s\n",
			        p.save_oxygen ? "on" : "off", p.save_iron ? "on" : "off");

			if (p.test_rt)
				fprintf(rf, "Mode:    %s test  tau_slab=%.4g\n",
				        p.test_mode, p.tau_slab);
			else if (strcmp(p.corona, "cutoffpl") == 0)
				fprintf(rf, "Corona:  cutoffpl  Gamma=%.4f  E_cut=%.2f keV\n",
				        p.Gamma, p.E_cut);
			else if (strcmp(p.corona, "nthcomp") == 0)
				fprintf(rf, "Corona:  nthcomp  Gamma=%.4f  kT_e=%.2f keV  kT_bb=%.4f keV\n",
				        p.Gamma, p.kT_e, p.kT_bb);
			else if (strcmp(p.corona, "comptt") == 0)
				fprintf(rf, "Corona:  comptt  kT_e=%.2f keV  kT_bb=%.4f keV  taup=%.2f\n",
				        p.kT_e, p.kT_bb, p.taup);
			else if (strcmp(p.corona, "blackbody") == 0)
				fprintf(rf, "Corona:  blackbody  kT_bb=%.4f keV\n", p.kT_bb);
			else
				fprintf(rf, "Corona:  %s  Gamma=%.4f\n", p.corona, p.Gamma);

			fprintf(rf, "Slab:    nh=%.4f  zeta=%.4f (xi=%.4e)  frac=%.4f  Afe=%.2f\n",
			        p.nh, p.zeta, xi, p.frac, p.Afe);
			fprintf(rf, "Derived: nH=%.4e  incident_corona_flux=%.4e erg/cm^2/s\n",
			        nH, Fx * (p.frac > 0 ? p.frac / (1.0 + p.frac) : 1.0));
			fprintf(rf, "Normalization: %s\n", "reflionx flux");
			fprintf(rf, "Depth grid: %s, %d cells\n", p.log_depth ? "log" : "tanh (k=3)",
			        RTGrids::ND_MID_DEFAULT);
			fprintf(rf, "Kernel:  %s\n",
			        p.angsca ? "angle-dependent (KernelCache)"
			                 : "angle-mean (avgKernelCache)");
			const bool failed = ferror(rf);
			if (fclose(rf) != 0 || failed) { fprintf(stderr, "Error writing RUN.txt\n"); exit(1); }
		}
	}

	return p;
}

void snap_incidence(ModelParams& par, const RTGrids& g)
{
	if (par.incidence == -2.0)
	{
		par.i_incidence = -1;
		dao_log::info("Incidence: isotropic over all downward angles (mu < 0)\n");
		return;
	}

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

	dao_log::info("Incidence:  cos(theta)=%.6f  (g.mu[%d])\n", par.incidence, par.i_incidence);
}
