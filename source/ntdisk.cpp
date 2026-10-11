#include "ntdisk.h"
#include "constants.h"

#include <cmath>
#include <stdexcept>

namespace ntdisk {
namespace {
// CGS constants; G and M_sun are explicit rather than hidden unit conversions.
constexpr long double G = 6.67430e-8L;
constexpr long double solar_mass_g = 1.98847e33L;
constexpr long double sigma_SB = 5.670374419e-5L;
constexpr long double pi = 3.1415926535897932384626433832795L;

void positive(double value, const char* name)
{
    if (!std::isfinite(value) || value <= 0)
        throw std::invalid_argument(name);
}

void validate_spin(double spin)
{
    if (!std::isfinite(spin) || std::abs(spin) >= 1)
        throw std::invalid_argument("NT spin must be finite and strictly between -1 and 1");
}

struct Factors { long double A, B, C, D; };

template<class Function>
long double refine(const Function& f, long double lo, long double hi,
                   long double fl, long double fm, long double fh,
                   long double whole, long double tolerance, int depth)
{
    const long double mid = (lo + hi) / 2;
    const long double f1 = f((lo + mid) / 2), f2 = f((mid + hi) / 2);
    const long double left = (mid - lo) * (fl + 4*f1 + fm) / 6;
    const long double right = (hi - mid) * (fm + 4*f2 + fh) / 6;
    const long double delta = left + right - whole;
    if (std::abs(delta) <= 15*tolerance)
        return left + right + delta/15;
    if (depth == 0)
        throw std::runtime_error("NT Eq. (8) quadrature did not converge");
    return refine(f, lo, mid, fl, f1, fm, left, tolerance/2, depth-1)
         + refine(f, mid, hi, fm, f2, fh, right, tolerance/2, depth-1);
}

Factors factors(double radius, double spin, double isco)
{
    const long double r = radius, a = spin;
    const long double u = 1/std::sqrt(r);
    // Hubeny & Hubeny (1998), Eqs. (5), (6), (7), respectively.
    Factors result{1-2*u*u+a*a*u*u*u*u,
                   1-3*u*u+2*a*u*u*u,
                   1-4*a*u*u*u+3*a*a*u*u*u*u, 0};
    if (radius == isco) return result;

    // Hubeny & Hubeny (1998), Eq. (8): D=(1/(2 sqrt(r))) int_ri^r ... dx.
    // Numerical change of variable ONLY: s=sqrt(x), y=log(s/s_i), y=L*t.
    // D=int_0^1 L*(s/sqrt(r))*(s^4-6s^2+8as-3a^2)/(s^4-3s^2+2as) dt.
    // Factor the numerator about s_i, where it vanishes by the ISCO condition
    // in Bardeen et al. (1972), Eq. (2.20). This avoids subtracting two nearly
    // equal integrals or large polynomial terms close to the inner edge.
    const long double si = std::sqrt(static_cast<long double>(isco));
    const long double ui = 1/si;
    const long double L = std::log1p((r-isco)/isco)/2;
    const auto integrand = [=](long double t) {
        const long double y = L*t;
        const long double v = std::exp(-y), d = -std::expm1(-y);
        const long double invs = ui*v;
        const long double numerator = d*((4-12*ui*ui+8*a*ui*ui*ui)*v*v*v
            + d*(6-6*ui*ui)*v*v + 4*d*d*v + d*d*d);
        const long double denominator = 1-3*invs*invs+2*a*invs*invs*invs;
        return L*std::exp(y-L)*numerator/denominator;
    };
    const long double f0 = integrand(0), fm = integrand(0.5L), f1 = integrand(1);
    const long double whole = (f0+4*fm+f1)/6;
    result.D = refine(integrand, 0, 1, f0, fm, f1, whole,
                      std::abs(whole)*1e-12L, 28);
    if (!(result.A > 0 && result.B > 0 && result.C > 0 && result.D > 0)
        || !std::isfinite(result.D))
        throw std::runtime_error("NT relativistic factors are nonpositive or nonfinite");
    return result;
}

struct Local { double isco; long double omega2; Factors correction; };

Local local(const SystemParameters& system, double radius)
{
    positive(system.mass_msun, "NT mass_msun must be positive and finite");
    positive(system.accretion_rate_g_s, "NT accretion_rate_g_s must be positive and finite");
    positive(radius, "NT radius_rg must be positive and finite");
    const double isco = isco_radius_rg(system.spin);
    if (radius < isco)
        throw std::invalid_argument("NT annulus must be at or outside ISCO");
    // Hubeny & Hubeny (1998), definitions after Eq. (8): R=r*GM/c^2.
    // omega2=GM/R^3 is a convenient prefactor, NOT the Kerr coordinate Omega^2.
    const long double gm = G*solar_mass_g*system.mass_msun;
    const long double R = radius*gm/(phys::c*phys::c);
    const long double omega2 = gm/R/R/R;
    if (!(omega2 > 0) || !std::isfinite(omega2))
        throw std::overflow_error("NT physical radius/frequency cannot be represented");
    return {isco, omega2, factors(radius, system.spin, isco)};
}

double finite_output(long double value)
{
    const double output = static_cast<double>(value);
    if (!std::isfinite(output) || output <= 0)
        throw std::overflow_error("NT positive output cannot be represented as a double");
    return output;
}
} // namespace

double isco_radius_rg(double spin)
{
    validate_spin(spin);
    if (spin == 0) return 6;
    const long double a = spin;
    // At tiny |a|, 3-Z1 in Eq. (2.21) loses precision. Evaluate the same root
    // using the polynomial ISCO condition in BPT Eq. (2.20), Newton iteration
    // from s=sqrt(6). This branch is a numerical evaluation, not new physics.
    if (std::abs(a) < 1e-4L) {
        long double s = std::sqrt(6.L);
        for (int i = 0; i < 5; ++i)
            s -= (s*s*s*s-6*s*s+8*a*s-3*a*a)/(4*s*s*s-12*s+8*a);
        return static_cast<double>(s*s);
    }
    // Bardeen, Press & Teukolsky (1972), ApJ 178, 347, Eq. (2.21).
    // URL: https://articles.adsabs.harvard.edu/pdf/1972ApJ...178..347B
    // Units GM/c^2. Positive a uses the direct-orbit minus sign; negative a
    // uses the retrograde plus sign, with disk angular momentum positive.
    const long double z1 = 1+std::cbrt(1-a*a)*(std::cbrt(1+a)+std::cbrt(1-a));
    const long double z2 = std::sqrt(3*a*a+z1*z1);
    return static_cast<double>(3+z2-std::copysign(
        std::sqrt((3-z1)*(3+z1+2*z2)), a));
}

AnnulusRadiation evaluate_radiation(const SystemParameters& system, double radius)
{
    const auto annulus = local(system, radius);
    if (radius == annulus.isco) return {radius, annulus.isco, 0, 0};
    // Hubeny & Hubeny (1998), ApJ 505, 558, Eq. (10).
    // URL for all HH98 equation numbers: https://arxiv.org/html/astro-ph/9804288
    // F=(3/(8*pi))*GM*dotM/R^3*(D/B) is ONE-face local intrinsic flux.
    const long double flux = 3/(8*pi)*system.accretion_rate_g_s*annulus.omega2
        *annulus.correction.D/annulus.correction.B;
    // Same Eq. (10): sigma_SB*T_eff^4=F. No color factor or redshift applied.
    return {radius, annulus.isco, finite_output(flux),
            finite_output(std::sqrt(std::sqrt(flux/sigma_SB)))};
}

SurfaceDensityDiagnostic evaluate_surface_density_diagnostic(
    const SystemParameters& system, double radius, const DiagnosticParameters& parameters)
{
    positive(parameters.alpha0, "NT diagnostic alpha0 must be positive and finite");
    positive(parameters.electron_opacity_cm2_g, "NT diagnostic electron opacity must be positive and finite");
    const auto annulus = local(system, radius);
    if (radius == annulus.isco)
        throw std::invalid_argument("NT surface-density diagnostic is singular at ISCO; use radius_rg>ISCO");
    const auto f = annulus.correction;
    const long double mdot = system.accretion_rate_g_s;
    const long double kappa = parameters.electron_opacity_cm2_g;
    const long double opacity_factor = kappa/(8*pi*phys::c);
    // HH98 Eq. (16), radiation-pressure approximation, with dotM^2 correction:
    // Pbar/rhobar=3*(GM/R^3)*dotM^2*(sigma_e/(8*pi*m_H*c))^2*D^2/(B*C).
    // Both arXiv v1 and the published ApJ Eq. (16) print dotM to the first power;
    // this is dimensionally inconsistent. The square follows from
    // Eqs. (14),(15),(18), or Pbar/rhobar
    // =Q*h^2/3 in Sec. 3.1 together with Eqs. (49),(53). It also reproduces
    // Eq. (17). This is an explicit algebraic correction, not a verbatim Eq. (16).
    // Replacing sigma_e/m_H by constant kappa is a stated composition
    // generalization. The default uses solar84 and n_e=1.21*n_H; explicitly
    // passing sigma_T/m_H instead reproduces the paper's pure-H choice.
    const long double pressure_density = 3*annulus.omega2*mdot*mdot
        *opacity_factor*opacity_factor*f.D*f.D/(f.B*f.C);
    // HH98 Eq. (18): m0=(16*pi/3)*(c/kappa)^2/(dotM*alpha0*sqrt(GM/R^3))
    //                   *B*C/(A*D).
    // alpha0 is the vertically averaged TOTAL-pressure prescription, Eq. (14).
    const long double m0 = 16*pi/3*(phys::c/kappa)*(phys::c/kappa)
        /(mdot*parameters.alpha0*std::sqrt(annulus.omega2))*f.B*f.C/(f.A*f.D);
    // HH98 text immediately after Eq. (12): m0=Sigma/2; BOTH-face diagnostic.
    // Thomson definition tau=int sigma_T*n_e dz. Using the fixed reference
    // n_e=1.21*n_H and rho=mass_per_H*n_H gives tau=kappa_T*column_mass.
    // HH98 Eq. (27) defines optical depth; the text after Eq. (12) gives
    // m0=Sigma/2. Hence the one-sided depth is kappa_T*m0, total is twice it.
    // Use reference kappa_T even if the analytic disk opacity is overridden.
    const long double tau_mid = reference_thomson_opacity_cm2_g*m0;
    return {finite_output(2*m0), finite_output(pressure_density),
            finite_output(tau_mid), finite_output(2*tau_mid)};
}
} // namespace ntdisk
