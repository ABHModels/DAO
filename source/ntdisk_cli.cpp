#include "ntdisk.h"

#include <cerrno>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void usage()
{
    std::cout << "Usage: ntdisk --mass-msun M --mdot-g-s MDOT --spin A --alpha0 ALPHA\n"
                 "              --radius-rg R [--radius-rg R ...] [--kappa-es CM2_PER_G]\n"
                 "NT zero-ISCO-stress local one-face radiation; radius must exceed ISCO.\n"
                 "Sigma is ONLY a radiation-pressure/electron-scattering diagnostic,\n"
                 "using vertically averaged stress = alpha0 * total pressure.\n"
                 "Default opacity: solar84 with n_e=1.21*n_H; mdot is g/s, not L/L_Edd.\n"
                 "Thomson depths use this reference composition: surface-to-midplane and total.\n";
}

double number(const char* text)
{
    errno = 0;
    char* end = nullptr;
    const double value = std::strtod(text, &end);
    if (end == text || *end != '\0' || errno == ERANGE)
        throw std::invalid_argument("invalid numerical argument");
    return value;
}
} // namespace

int main(int argc, char** argv)
{
    try {
        ntdisk::SystemParameters system{};
        ntdisk::DiagnosticParameters diagnostic{};
        bool mass = false, mdot = false, spin = false, alpha = false;
        std::vector<double> radii;
        for (int i = 1; i < argc; ++i) {
            const std::string option = argv[i];
            if (option == "--help") { usage(); return 0; }
            if (i+1 == argc) throw std::invalid_argument("missing option value");
            const double value = number(argv[++i]);
            if (option == "--mass-msun") { system.mass_msun = value; mass = true; }
            else if (option == "--mdot-g-s") { system.accretion_rate_g_s = value; mdot = true; }
            else if (option == "--spin") { system.spin = value; spin = true; }
            else if (option == "--alpha0") { diagnostic.alpha0 = value; alpha = true; }
            else if (option == "--radius-rg") radii.push_back(value);
            else if (option == "--kappa-es") diagnostic.electron_opacity_cm2_g = value;
            else throw std::invalid_argument("unknown option: " + option);
        }
        if (!mass || !mdot || !spin || !alpha || radii.empty())
            throw std::invalid_argument("mass, mdot, spin, alpha0 and radius are required");
        struct Row {
            ntdisk::AnnulusRadiation radiation;
            ntdisk::SurfaceDensityDiagnostic diagnostic;
        };
        std::vector<Row> rows;
        for (double radius : radii)
            rows.push_back({ntdisk::evaluate_radiation(system, radius),
                ntdisk::evaluate_surface_density_diagnostic(system, radius, diagnostic)});
        std::cout << "radius_rg,isco_rg,Teff_K,F_intrinsic_erg_cm2_s,Sigma_diagnostic_g_cm2,"
                     "Pbar_over_rhobar_diagnostic_cm2_s2,tau_Thomson_mid_diagnostic,"
                     "tau_Thomson_total_diagnostic\n";
        std::cout << std::scientific << std::setprecision(16);
        for (const auto& row : rows)
            std::cout << row.radiation.radius_rg << ',' << row.radiation.isco_rg << ','
                      << row.radiation.effective_temperature_K << ','
                      << row.radiation.intrinsic_flux_erg_cm2_s << ','
                      << row.diagnostic.surface_density_g_cm2 << ','
                      << row.diagnostic.mean_pressure_over_density_cm2_s2 << ','
                      << row.diagnostic.thomson_depth_midplane << ','
                      << row.diagnostic.thomson_depth_total << '\n';
    } catch (const std::exception& error) {
        std::cerr << "ntdisk: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
