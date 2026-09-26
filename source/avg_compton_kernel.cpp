#include "avg_compton_kernel.h"
#include "compton_kernel.h"          // profil_exact / profil_exact_ap
#include "compton_cross_section.h"
#include "kernel_row_spool.h"
#include "rt_parallel.h"
#include "constants.h"
#include "kernel_quadrature.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

// ============================================================
// Gauss-Legendre quadrature over the scattering-angle cosine
//
//   K̄(x; x₁) = 2π ∫₋₁¹ R(x, x₁, cosθ) d(cosθ)
//
// Adaptive quadrature resolves the narrow low-temperature recoil peak.
// ============================================================

// ============================================================
// Angle-mean kernel element
//
//   K̄(x; x₁, T) = 2π ∫₋₁¹ R(x, x₁, cosθ; 1/Θ) d(cosθ)
//
// ktype == 0 → approximate redistribution (profil_exact_ap),
// otherwise  → exact redistribution (profil_exact).
// ============================================================
double angle_mean_kernel_element(double x, double x1, double T_K,
                                 const int ktype)
{
	const double boltz_eV = 8.617333262e-5;   // Boltzmann constant [eV/K]
	const double mec2_eV  = 511.0e3;          // electron rest energy [eV]

	double x_inv = mec2_eV / (boltz_eV * T_K);   // m_e c² / kT
	const double delta=x-x1;
	const double peak_q=fabs(delta);
	const double peak_cos=1.0-peak_q/(x*x1);
	if (ktype!=0) {
		const double q=std::min(2.0*x*x1,peak_q);
		const double gamma=q<=0 ? 1.0 :
			0.5*(delta+sqrt(delta*delta+2*q)*sqrt(1.0+2.0/q));
		if ((gamma-1.0)*x_inv>120.0) return 0.0;
	}

	auto profile=[&](double costh) {
		const double c=std::min(std::nextafter(1.0,0.0),
		                        std::max(-1.0,costh));
		return ktype==0 ? profil_exact_ap(x,x1,c,x_inv)
		                : profil_exact(x,x1,c,x_inv);
	};
	std::array<double,5> cuts{-1.0,1.0};
	size_t count=2;
	if (peak_cos>=-1.0 && peak_cos<=1.0) {
		const double dc=4.0*sqrt(2.0/x_inv)/std::max(x,x1);
		cuts[count++]=std::max(-1.0,peak_cos-dc);
		cuts[count++]=peak_cos;
		cuts[count++]=std::min(1.0,peak_cos+dc);
	}
	const double thermal=x*sqrt(2.0/x_inv);
	const double integral=fabs(x1-x)<4.0*thermal
		? kernel_quadrature::integrate_fixed<96>(profile,cuts,count)
		: kernel_quadrature::integrate_fixed<48>(profile,cuts,count);
	return 2.0*M_PI*integral;
}

static double angle_mean_energy_cell(double x,double lo,double hi,double T_K,int ktype)
{
	auto profile=[&](double x1) {
		return angle_mean_kernel_element(x,x1,T_K,ktype);
	};
	const double thermal=x*sqrt(2.0*8.617333262e-5*T_K/511000.0);
	std::array<double,6> cuts{lo,hi};
	size_t count=2;
	auto add=[&](double point) {
		if(point>lo && point<hi) cuts[count++]=point;
	};
	add(x);
	add(x-4*thermal);
	add(x+4*thermal);
	const double denominator=1.0/x-2.0;
	if(denominator>0) add(1.0/denominator);
	const double integral=lo<=x && x<=hi
		? kernel_quadrature::integrate_fixed<128>(profile,cuts,count)
		: kernel_quadrature::integrate_fixed<48>(profile,cuts,count);
	return integral/(hi-lo);
}

// ============================================================
// avgKernelCache implementation — banded storage, detailed balance
// (upper triangle ne1 >= ne only)
//
// Binary cache file format (version 07):
//   magic       [8 bytes]  "AVKRN07\0"
//   NT, NE      [2×4 bytes]
//   data_size   [8 bytes]
//   T_grid      [NT doubles]
//   x_grid      [NE doubles]
//   theta       [NT doubles]
//   band_lo     [NT*NE ints]
//   band_hi     [NT*NE ints]
//   band_off    [NT*NE longs]
//   glo         [NT*NE ints]
//   ghi         [NT*NE ints]
//   data        [data_size doubles]
// ============================================================
static const char CACHE_MAGIC[8] = "AVKRN07";

void avgKernelCache::write_header(FILE* fp) const
{
	long n_ge = long(NT) * NE;

	fwrite(CACHE_MAGIC, 1, 8, fp);
	fwrite(&NT,        sizeof(int),  1, fp);
	fwrite(&NE,        sizeof(int),  1, fp);
	fwrite(&data_size, sizeof(long), 1, fp);
	fwrite(T_grid,   sizeof(double), NT,        fp);
	fwrite(x_grid,   sizeof(double), NE,        fp);
	fwrite(theta,    sizeof(double), NT,        fp);
	fwrite(band_lo,  sizeof(int),    n_ge,      fp);
	fwrite(band_hi,  sizeof(int),    n_ge,      fp);
	fwrite(band_off, sizeof(long),   n_ge,      fp);
	fwrite(glo,      sizeof(int),    n_ge,      fp);
	fwrite(ghi,      sizeof(int),    n_ge,      fp);
}

void avgKernelCache::save(const char* filename) const
{
	KernelCacheOutput output(filename);
	FILE* fp = output.file();
	write_header(fp);
	fwrite(data,     sizeof(double), data_size, fp);
	output.commit();

	fprintf(stdout, "  avgKernelCache: saved to %s (%.1f MB)\n",
	        filename, data_size * 8.0 / (1024.0 * 1024.0));
}

bool avgKernelCache::load(const char* filename, const double* expected_ene_eV)
{
	// NT, NE, T_grid must be set before calling.
	FILE* fp = fopen(filename, "rb");
	if (!fp) return false;

	char magic[8];
	if (fread(magic, 1, 8, fp) != 8 || memcmp(magic, CACHE_MAGIC, 8) != 0)
	{
		fprintf(stderr,
			"  WARNING: cache file %s has wrong magic (expected \"%s\", got \"%.7s\")\n"
			"           File may be from an older format version. Will recompute.\n",
			filename, CACHE_MAGIC, magic);
		fclose(fp); return false;
	}

	int fNT, fNE;
	fread(&fNT, sizeof(int), 1, fp);
	fread(&fNE, sizeof(int), 1, fp);
	if (fNT != NT || fNE != NE)
	{
		fprintf(stderr,
			"  WARNING: cache file %s has mismatched dimensions:\n"
			"           file: NT=%d NE=%d\n"
			"           need: NT=%d NE=%d\n"
			"           Will recompute.\n",
			filename, fNT, fNE, NT, NE);
		fclose(fp); return false;
	}

	long f_data_size;
	fread(&f_data_size, sizeof(long), 1, fp);

	// Verify T_grid
	double* fT = new double[NT];
	fread(fT, sizeof(double), NT, fp);
	bool match = true;
	for (int i = 0; i < NT; ++i)
		if (fabs(fT[i] - T_grid[i]) > 1e-6 * fT[i]) { match = false; break; }
	delete[] fT;
	if (!match)
	{
		fprintf(stderr,
			"  WARNING: cache file %s has mismatched T_grid. Will recompute.\n",
			filename);
		fclose(fp); return false;
	}

	// Read x_grid and theta
	x_grid = new double[NE];
	fread(x_grid, sizeof(double), NE, fp);
	if (expected_ene_eV) {
		bool grid_match=true;
		for (int i=0; i<NE; ++i)
			if (fabs(x_grid[i]-expected_ene_eV[i]/511000.0) >
			    1e-12*expected_ene_eV[i]/511000.0) { grid_match=false; break; }
		if (!grid_match) {
			delete[] x_grid; x_grid=nullptr;
			fclose(fp); return false;
		}
	}
	theta = new double[NT];
	fread(theta, sizeof(double), NT, fp);

	long n_ge = long(NT) * NE;

	band_lo  = new int[n_ge];
	band_hi  = new int[n_ge];
	band_off = new long[n_ge];
	fread(band_lo,  sizeof(int),  n_ge, fp);
	fread(band_hi,  sizeof(int),  n_ge, fp);
	fread(band_off, sizeof(long), n_ge, fp);

	glo = new int[n_ge];
	ghi = new int[n_ge];
	fread(glo, sizeof(int), n_ge, fp);
	fread(ghi, sizeof(int), n_ge, fp);

	data_size = f_data_size;
	data = payload.map(fp,size_t(data_size));
	if(data) {
		fclose(fp);
		fprintf(stdout,"  avgKernelCache: mapped %s (%.1f MB virtual payload, demand-paged)\n",
		        filename,data_size*8.0/(1024.0*1024.0));
		return true;
	}
	data = new double[data_size];
	size_t nread = fread(data, sizeof(double), data_size, fp);
	fclose(fp);

	if ((long)nread != data_size) return false;

	double fill = 100.0 * data_size / (double(NT) * NE * NE);
	fprintf(stdout, "  avgKernelCache: loaded from %s (%.1f MB, %.1f%% dense fill)\n",
	        filename, data_size * 8.0 / (1024.0 * 1024.0), fill);
	return true;
}

void avgKernelCache::init(int n_ene, const double* ene_eV,
                          int /*n_ang*/, const double* /*mu*/, const double* /*wt*/,
                          const int ktype,
                          int nt_user, const double* T_user)
{
	NE = n_ene;

	// Temperature grid: caller-supplied (e.g. single isothermal value in
	// test mode) or the default log-spaced N_T_CACHE grid for production.
	if (nt_user > 0 && T_user)
	{
		NT     = nt_user;
		T_grid = new double[NT];
		for (int i = 0; i < NT; ++i)
			T_grid[i] = T_user[i];
	}
	else
	{
		NT     = N_T_CACHE;
		T_grid = new double[NT];
		double log_lo = log10(T_CACHE_LO);
		double log_hi = log10(T_CACHE_HI);
		double dlog   = (NT > 1) ? (log_hi - log_lo) / (NT - 1) : 0.0;
		for (int i = 0; i < NT; ++i)
			T_grid[i] = pow(10.0, log_lo + i * dlog);
	}

	// Cache filenames
	const char* cache_dir = getenv("COMPTON_CACHE_DIR");
	if (!cache_dir) cache_dir = ".";

	const char* ksuffix = (ktype == 0) ? "_ap" : "";
	char norm_file[512], raw_file[512];
	snprintf(norm_file, sizeof(norm_file),
	         "%s/avgkernel_norm_NE%d_NT%d%s.bin",
	         cache_dir, NE, NT, ksuffix);
	snprintf(raw_file, sizeof(raw_file),
	         "%s/avgkernel_NE%d_NT%d%s.bin",
	         cache_dir, NE, NT, ksuffix);

	if (load(norm_file,ene_eV))
	{
		fprintf(stderr, "\n Load normalized angle-mean Compton kernel");
		return;
	}
	if (load(raw_file,ene_eV))
	{
		fprintf(stderr,
			"\n  WARNING: Loaded un-normalized angle-mean kernel %s\n"
			"           Normalized kernel %s not found.\n",
			raw_file, norm_file);
		return;
	}

	fprintf(stdout, "  avgKernelCache: no valid cache, computing "
	        "(banded, upper-triangle, detailed balance)...\n");

	// Convert energy grid to dimensionless x = E/(m_e c²)
	const double mec2_eV  = 511.0e3;
	const double boltz_eV = 8.617333262e-5;
	x_grid = new double[NE];
	for (int ie = 0; ie < NE; ++ie)
		x_grid[ie] = ene_eV[ie] / mec2_eV;

	// Dimensionless temperature θ = kT/(m_e c²)
	theta = new double[NT];
	for (int i = 0; i < NT; ++i)
		theta[i] = boltz_eV * T_grid[i] / mec2_eV;

	// Allocate band structure
	long n_ge = long(NT) * NE;
	band_lo  = new int[n_ge];
	band_hi  = new int[n_ge];
	band_off = new long[n_ge];

	KernelRowSpool computed_rows;

	const double KMIN_ABS = 1e-30;
	const double KMIN_REL = 1e-10;

	// --- Pass 1: find band limits (upper triangle only: ne1 >= ne) ---
	fprintf(stdout, "  avgKernelCache: pass 1 — finding band limits (upper triangle)...\n");

	computed_rows.evaluate(size_t(n_ge),size_t(NE),[&](size_t index,double* row_buf,int& lo_ne1,int& hi_ne1) {
		const long r=long(index);
		const int ne=int(index%size_t(NE));
		const int iT=int(index/size_t(NE));
		lo_ne1=NE; hi_ne1=-1;
			double K_max = 0.0;

			// Only scan ne1 >= ne (upper triangle)
			for (int ne1 = ne; ne1 < NE; ++ne1)
			{
				double val;
				const double lo=NE<2 ? x_grid[ne1] :
					(ne1 ? 0.5*(x_grid[ne1-1]+x_grid[ne1])
					     : x_grid[ne1]-0.5*(x_grid[ne1+1]-x_grid[ne1]));
				const double hi=NE<2 ? x_grid[ne1] :
					(ne1+1<NE ? 0.5*(x_grid[ne1]+x_grid[ne1+1])
					          : x_grid[ne1]+0.5*(x_grid[ne1]-x_grid[ne1-1]));
				const double thermal=x_grid[ne1]*sqrt(2.0*theta[iT]);
				const double edge_margin=8.0*thermal;
				// Average near the forward/back recoil edges, where point
				// samples alias a narrow profile on the transfer grid.
				const double recoil_den=1.0/x_grid[ne]-2.0;
				const bool near_forward=lo<=x_grid[ne]+edge_margin &&
				                        hi>=x_grid[ne]-edge_margin;
				const bool near_back=recoil_den>0 &&
					lo<=1.0/recoil_den+edge_margin && hi>=1.0/recoil_den-edge_margin;
				if(NE>1 && hi-lo>2*thermal && (near_forward || near_back))
					val=angle_mean_energy_cell(x_grid[ne],std::max(1e-30,lo),hi,
				                              T_grid[iT],ktype);
				else
					val=angle_mean_kernel_element(x_grid[ne],x_grid[ne1],T_grid[iT],ktype);
				row_buf[ne1] = val;
				if (val > KMIN_ABS)
				{
					if (ne1 < lo_ne1) lo_ne1 = ne1;
					hi_ne1 = ne1;
				}
				if (val > K_max) K_max = val;
			}

			if (hi_ne1 < 0) { lo_ne1 = 0; hi_ne1 = -1; }

			// Trim band edges with dynamic threshold
			if (hi_ne1 >= 0)
			{
				double thresh = K_max * KMIN_REL;

				while (lo_ne1 <= hi_ne1)
				{
					if (row_buf[lo_ne1] > thresh) break;
					++lo_ne1;
				}
				while (hi_ne1 >= lo_ne1)
				{
					if (row_buf[hi_ne1] > thresh) break;
					--hi_ne1;
				}
				if (hi_ne1 < lo_ne1) { lo_ne1 = 0; hi_ne1 = -1; }
			}

			band_lo[r] = lo_ne1;
			band_hi[r] = hi_ne1;
	});

	// Compute per-(iT,ne) global band including both triangles.
	// Upper triangle: directly from stored band.
	// Lower triangle: ne1 < ne where K̄(ne1, ne) is stored,
	//   i.e., ne is within the band of row (iT, ne1).
	glo = new int[n_ge];
	ghi = new int[n_ge];
	for (int iT = 0; iT < NT; ++iT)
	for (int ne = 0; ne < NE; ++ne)
	{
		long r = row(iT, ne);
		int lo_g = (band_hi[r] >= band_lo[r]) ? band_lo[r] : NE;
		int hi_g = band_hi[r];

		// Lower triangle (via detailed balance): scan ne1 < ne downward
		for (int ne1 = ne - 1; ne1 >= 0; --ne1)
		{
			long rr = row(iT, ne1);
			if (ne >= band_lo[rr] && ne <= band_hi[rr])
			{
				lo_g = ne1;
				if (hi_g < ne1) hi_g = ne1;
			}
		}

		glo[r] = (hi_g >= 0 || lo_g < NE) ? lo_g : 0;
		ghi[r] = hi_g;
	}

	// Compute offsets and total data size
	data_size = 0;
	for (long r = 0; r < n_ge; ++r)
	{
		band_off[r] = data_size;
		int bw = band_hi[r] - band_lo[r] + 1;
		if (bw > 0) data_size += long(bw);
	}

	double fill = 100.0 * data_size / (double(NT) * NE * NE);
	fprintf(stdout, "  avgKernelCache: data_size = %ld doubles (%.1f MB, %.1f%% dense fill)\n",
	        data_size, data_size * 8.0 / (1024.0 * 1024.0), fill);

	// Reuse the exact values evaluated during the band scan.
	fprintf(stdout, "  avgKernelCache: normalizing retained rows one temperature at a time...\n");
	KernelCacheOutput output(norm_file);
	write_header(output.file());
	const off_t payload_offset=::ftello(output.file());


	// --- Normalize via the angle-mean A23 sum rule ---
	// For each incoming energy x₁, the kernel integrated over the
	// outgoing energy must reproduce the exact Compton cross section:
	//
	//   σ(x₁, T)/σ_T = (1/x₁) ∫ dx x · K̄(x; x₁)
	//
	// The K() accessor spans both triangles (direct + detailed balance),
	// so the integral covers the full energy range correctly.
	fprintf(stdout, "  Normalizing angle-mean kernel via A23 sum rule...\n");

	for (int iT = 0; iT < NT; ++iT)
	{
		fprintf(stdout, "    T=%.2e K (%d/%d)...", T_grid[iT], iT + 1, NT);
		fflush(stdout);
		data=computed_rows.view(payload,size_t(data_size));

		std::vector<double> norm(NE);
		for (int iteration=0; iteration<20; ++iteration) {
		rt_parallel_depths(NE, [&](int nin) { // independent incoming energies
			double x1 = x_grid[nin];

			double integral = 0.0;
			for (int nout = 1; nout < NE; ++nout)   // nout = outgoing energy x
			{
				double dx  = x_grid[nout] - x_grid[nout - 1];
				double f_a = x_grid[nout - 1] * K(iT, nout - 1, nin);
				double f_b = x_grid[nout]     * K(iT, nout,     nin);
				integral += 0.5 * (f_a + f_b) * dx;
			}

			double sigma_raw   = integral / x1;
			double sigma_exact = compton_cross_section(ene_eV[nin], T_grid[iT])
			                   / phys::sigma_T;

				norm[nin] = (sigma_raw > 1e-30) ? sigma_exact / sigma_raw : 0.0;
		}, "DAO_KERNEL_THREADS", 8);
		double max_error=0;
		for (double ratio:norm) {
			if (!(ratio>0) || !std::isfinite(ratio))
				throw std::runtime_error("kernel: zero or nonfinite scattering normalization");
			max_error=std::max(max_error,fabs(ratio-1.0));
		}
		if (max_error<1e-3) break;

		// Symmetric scaling preserves detailed balance on both triangles.
		for (int ne_out = 0; ne_out < NE; ++ne_out)
		{
			long r  = row(iT, ne_out);
			int  lo = band_lo[r];
			int  hi = band_hi[r];
			if (hi < lo) continue;

			double* dst = data + band_off[r];
			for (int ne_in = lo; ne_in <= hi; ++ne_in)
				dst[ne_in - lo] *= sqrt(norm[ne_out]*norm[ne_in]);
		}
		if (iteration==19)
			throw std::runtime_error("kernel: scattering normalization did not converge");
		}
		const long begin=band_off[long(iT)*NE];
		const long end=(iT+1<NT)?band_off[long(iT+1)*NE]:data_size;
		if(end>begin && fwrite(data+begin,sizeof(double),size_t(end-begin),output.file())!=size_t(end-begin))
			throw std::runtime_error("kernel: normalized temperature write failed");
		fprintf(stdout, " done.\n");
	}

	// Save normalized kernel
	data=output.map_payload(payload,size_t(data_size),payload_offset);
	output.commit();
	fprintf(stdout, "  Saved normalized angle-mean kernel to %s\n", norm_file);
}

void avgKernelCache::free_memory()
{
	delete[] T_grid;   T_grid   = nullptr;
	delete[] x_grid;   x_grid   = nullptr;
	delete[] theta;    theta    = nullptr;
	delete[] band_lo;  band_lo  = nullptr;
	delete[] band_hi;  band_hi  = nullptr;
	delete[] band_off; band_off = nullptr;
	delete[] glo;      glo      = nullptr;
	delete[] ghi;      ghi      = nullptr;
	payload.release(data);
	NT = NE = 0;
	data_size = 0;
}

int avgKernelCache::find_T(double T_K) const
{
	double logT = log10(T_K);
	int lo = 0, hi = NT - 1;
	while (lo < hi - 1)
	{
		int mid = (lo + hi) / 2;
		if (log10(T_grid[mid]) < logT) lo = mid;
		else hi = mid;
	}
	if (fabs(log10(T_grid[lo]) - logT) <= fabs(log10(T_grid[hi]) - logT))
		return lo;
	return hi;
}
