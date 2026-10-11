#ifndef DAO_THERMAL_BALANCE_H
#define DAO_THERMAL_BALANCE_H

#include "avg_compton_kernel.h"
#include "compton_kernel.h"
#include <vector>
#include <map>
#include <memory>

namespace dao_thermal {

class DirectionalScattering;

// Scattering operator at the actual gas temperature, in either angular mode.
// Interpolate the retained upper triangle in log T, reconstruct the reverse
// transitions at T, then symmetrically normalize to the same opacity used by RT.
// Pair scaling preserves detailed balance; separate row normalization would not.
// electron_density is the local free-electron density [cm^-3] from Cloudy.
// Opacity, redistribution and energy exchange all use this same density.
class Scattering {
    struct Pair { int i,j; double down,up; };
    std::vector<Pair> pairs;
    std::vector<double> opacity_, exchange_;
    std::shared_ptr<const DirectionalScattering> directional_;
public:
    Scattering() = default;
    Scattering(const avgKernelCache&, double temperature, double electron_density);
    Scattering(const KernelCache&, double temperature, double electron_density,
               int angles, const double* angular_weights);
    const std::vector<double>& opacity() const { return opacity_; }
    int angles() const;
    // Field/emission arrays are NE for angle-mean, or [angle][NE] for directional.
    void emissivity(const double* J, double* j) const;
    double energy_gain(const double* J) const;
    double normalization_error = 0;
};

std::vector<double> energy_weights(int n, const double* energy);

struct Budget {
    double absorption=0, emission=0, compton=0;
    double residual() const { return absorption-emission+compton; }
    double scale() const;
    double relative() const;
};
Budget budget(int n, const double* weights, const double* J,
              const double* absorption, const double* emission,
              const Scattering& scattering, const double* intensity=nullptr);

// Response of one constant cell to fixed incoming rays from its neighbours.
// Scattering is solved within the cell in the same angular mode as global RT.
// The stable residual uses k*phi*(I_in-S), not subtraction of trapped terms.
struct CellResponse {
    std::vector<double> mean, intensity;
    double residual=0, scale=0;
    double relative() const;
};
CellResponse local_response(int ne, int na, const double* energy_weights,
                            const double* mu, const double* angular_weights,
                            double width, const double* absorption,
                            const double* emission, const Scattering& scattering,
                            const double* incoming, const double* seed,
                            int max_iterations=500);

// Brackets a stable root (positive heating below, negative heating above) in
// log T. All samples must come from ONE fixed radiation/column context (either
// fixed J in a diagnostic, or fixed incoming cell rays in production).
// Exhausting the table or the evaluation budget is an error, never convergence.
class TemperatureRoot {
    double initial_, minimum_, maximum_;
    std::map<double,double> samples_;
    int expansion_=0;
public:
    TemperatureRoot(double initial, double minimum, double maximum);
    void record(double temperature, double residual);
    double next();
};

} // namespace dao_thermal
#endif
