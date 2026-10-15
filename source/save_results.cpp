#include "run_log.h"
#include "save_results.h"
#include "rt_grids.h"
#include "cddefines.h"
#include "cddrive.h"
#include "lines.h"
#include "rfield.h"
#include "dense.h"
#include "iso.h"
#include "taulines.h"
#include "transition.h"
#include "yield.h"
#include "heavy.h"
#include "ionbal.h"
#include "physconst.h"
#include <cstdio>
#include <cmath>
#include <string>
#include <sys/stat.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>

// Cloudy redefines fopen; undo it for file output
#undef fopen

namespace {
std::filesystem::path ion_path(const ModelParams& par, const char* folder)
{
    return std::filesystem::path(par.run_dir[0] ? par.run_dir : "results") / folder / "ion_fractions.dat";
}

template<size_t N> void write_ions(const RadField& rad, const RTGrids& g,
                                  const ModelParams& par, int iter, const char* folder,
                                  const char* symbol, const std::vector<IonFractions<N>>& ions)
{
    if(ions.size()!=size_t(g.ND_MID)) throw std::runtime_error("Missing final Cloudy ion fractions");
    for(const auto& row:ions) {
        if(row.iteration!=iter) throw std::runtime_error("Stale Cloudy ion fractions");
        if(!std::isfinite(row.gas_density) || row.gas_density<0)
            throw std::runtime_error("Invalid final elemental density");
        for(double value:row.fraction)
            if(!std::isfinite(value) || value<0 || value>1+1e-5)
                throw std::runtime_error("Invalid final ion fraction");
    }
    const auto path=ion_path(par,folder);
    const std::string temporary=path.string()+".tmp";
    std::filesystem::create_directories(path.parent_path());
    try {
        std::ofstream out;
        out.exceptions(std::ios::failbit|std::ios::badbit);
        out.open(temporary);
        out<<std::scientific<<std::setprecision(16)
           <<"# "<<symbol<<" ion fractions at final convergence; hash="<<par.run_hash<<" iteration="<<iter<<'\n'
           <<"# Atomic state from the last accepted-temperature Cloudy pass (preceding outer J); no diagnostic rerun.\n"
           <<"# f(X^q+) = n(X^q+)/n_X,gas, linear fractions; no renormalization over ion stages.\n"
           <<"# Molecules can make the sum below 1; absent elements have density=0 and all fractions=0.\n"
           <<"# depth_index is zero-based; z_cm is the DAO cell midpoint measured inward from the top.\n"
           <<"# tau_ref uses n_e,ref=1.21*nH; n_e_cm-3 is the local Cloudy electron density.\n"
           <<"# Ion columns run from neutral (I, q=0) to fully stripped; spectroscopic stage = q+1.\n"
           <<"# depth_index z_cm tau_ref T_K n_e_cm-3 n_"<<symbol<<"_gas_cm-3";
        const char* roman[]={"I","II","III","IV","V","VI","VII","VIII","IX",
            "X","XI","XII","XIII","XIV","XV","XVI","XVII","XVIII","XIX",
            "XX","XXI","XXII","XXIII","XXIV","XXV","XXVI","XXVII"};
        for(size_t q=0;q<N;++q) out<<' '<<symbol<<'_'<<roman[q];
        out<<'\n';
        double z=0;
        for(int d=0;d<g.ND_MID;++d) {
            out<<d<<' '<<z+0.5*g.dr[d]<<' '<<g.tau_mid[d]<<' '<<rad.T_K[d]
               <<' '<<rad.n_e[d]<<' '<<ions[d].gas_density;
            for(double value:ions[d].fraction) out<<' '<<value;
            out<<'\n'; z+=g.dr[d];
        }
        out.close();
        std::filesystem::rename(temporary,path);
    } catch(...) {
        std::error_code ignored;
        std::filesystem::remove(temporary,ignored);
        throw;
    }
    dao_log::info("Ion fractions: %s\n",path.c_str());
}
}

void clear_ion_fraction_output(const ModelParams& par)
{
    // Identical physics reuses the same hash, even when these output flags differ.
    // Remove only generated diagnostics so a failed rerun cannot expose old ions.
    for(const char* folder:{"O","Iron"}) {
        const auto path=ion_path(par,folder);
        std::filesystem::remove(path);
        std::filesystem::remove(path.string()+".tmp");
    }
}

void save_ion_fractions(const RadField& rad, const RTGrids& g,
                        const ModelParams& par, int iter)
{
    if(par.save_oxygen) write_ions(rad,g,par,iter,"O","O",rad.oxygen);
    if(par.save_iron) write_ions(rad,g,par,iter,"Iron","Fe",rad.iron);
}

void save_results(const RadField& rad, const RTGrids& g,
                  const ModelParams& par, int iter)
{
	const char* dir = par.run_dir[0] ? par.run_dir : "results";

	mkdir("results", 0755);
	if (par.run_dir[0])
		mkdir(dir, 0755);

	char fname[256];

	// 1. Emergent angle-dependent intensity at the actual upper slab face.
	{
		snprintf(fname, sizeof(fname), "%s/emergent_iter%03d.dat", dir, iter);
		FILE* fp = fopen(fname, "w");
		fprintf(fp, "# Emergent specific intensity at upper slab face  iter=%d\n", iter);
		fprintf(fp, "# Col 1: E [eV]\n");
		fprintf(fp, "# Col 2: F_corona (incident normal flux, stored as I_corona) [erg cm^-2 s^-1 eV^-1]\n");
		fprintf(fp, "# Col 3: I_disk (incident) [erg cm^-2 s^-1 eV^-1 sr^-1]\n");
		for (int nm = 0; nm < g.NA; ++nm)
			fprintf(fp, "# Col %d: I(mu=%.6f) [erg cm^-2 s^-1 eV^-1 sr^-1]\n",
			        nm + 4, g.mu[nm]);
		for (int ie = 0; ie < g.NE; ++ie)
		{
			fprintf(fp, "%.6e  %.6e  %.6e",
			        g.ene[ie], rad.illum.I_corona[ie], rad.illum.I_disk[ie]);
			for (int nm = 0; nm < g.NA; ++nm)
				fprintf(fp, "  %.6e", rad.Inu_top[nm][ie]);
			fprintf(fp, "\n");
		}
		fclose(fp);
		dao_log::detail("  Saved: %s\n", fname);
	}

	// 2. Angular moments at all depths
	{
		snprintf(fname, sizeof(fname), "%s/moments_iter%03d.dat", dir, iter);
		FILE* fp = fopen(fname, "w");
		if(!par.test_rt)
			fprintf(fp, "# Cell-volume-average angular moments: J0, J2, J3 (RT thermal balance) iter=%d\n", iter);
		else
			fprintf(fp, "# Cell-center angular moments: J0, J2, J3  iter=%d\n", iter);
		fprintf(fp, "# Col 1: E [eV]\n");
		fprintf(fp, "# Col 2: depth index\n");
		fprintf(fp, "# Col 3: tau_mid [reference Thomson depth, n_e,ref = 1.21*nH]\n");
		fprintf(fp, "# Col 4: T_K [K]\n");
		fprintf(fp, "# Col 5: J0 [erg cm^-2 s^-1 eV^-1]\n");
		fprintf(fp, "# Col 6: J2\n");
		fprintf(fp, "# Col 7: J3\n");
		for (int id = 0; id < g.ND_MID; ++id)
		for (int ie = 0; ie < g.NE; ++ie)
		{
			fprintf(fp, "%.6e  %d  %.6e  %.6e  %.6e  %.6e  %.6e\n",
			        g.ene[ie], id, g.tau_mid[id], rad.T_K[id],
			        rad.J0[id][ie], rad.J2[id][ie], rad.J3[id][ie]);
		}
		fclose(fp);
		dao_log::detail("  Saved: %s\n", fname);
	}

	// 3. Temperature and opacity profile
	{
		snprintf(fname, sizeof(fname), "%s/profile_iter%03d.dat", dir, iter);
		FILE* fp = fopen(fname, "w");
		fprintf(fp, "# Depth profile  iter=%d\n", iter);
		if(!par.test_rt)
			fprintf(fp, "# Fixed-T Cloudy: heating/cooling are diagnostics, not the temperature residual; see thermal_budget files.\n");
		fprintf(fp, "# Col 1: depth index\n");
		fprintf(fp, "# Col 2: tau_mid [reference Thomson depth, n_e,ref = 1.21*nH]\n");
		fprintf(fp, "# Col 3: T_K [K]\n");
		fprintf(fp, "# Col 4: n_e [cm^-3]\n");
		fprintf(fp, "# Col 5: Cloudy heating diagnostic [erg cm^-3 s^-1]\n");
		fprintf(fp, "# Col 6: cooling [erg cm^-3 s^-1]\n");
		fprintf(fp, "# Col 7: log_xi = log10(4piJ) [erg cm^-2 s^-1] / nh\n");
		fprintf(fp, "# Col 8: line_heat_Pdest [erg cm^-3 s^-1]\n");
		fprintf(fp, "# Col 9: heating_plus_new_line_heat = heating + newly computed line_heat_Pdest [erg cm^-3 s^-1]\n");
		for (int id = 0; id < g.ND_MID; ++id)
		{
			fprintf(fp, "%d  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e\n",
			        id, g.tau_mid[id], rad.T_K[id], rad.n_e[id],
			        rad.heating[id], rad.cooling[id], rad.log_xi[id],
			        rad.line_heat[id], rad.heating[id] + rad.line_heat[id]);
		}
		fclose(fp);
		dao_log::detail("  Saved: %s\n", fname);
	}

}
