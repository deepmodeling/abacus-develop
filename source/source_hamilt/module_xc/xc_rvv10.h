#ifndef ABACUS_XC_RVV10_H
#define ABACUS_XC_RVV10_H

#include "source_base/cubic_spline.h"

#include <array>

namespace Rvv10
{
constexpr int channel_count = 20;
constexpr double density_cutoff = 1.e-12;
constexpr double q_min = 1.e-4;
constexpr double q_cut = 0.5;
constexpr double rvv10_b_default = 6.3;
constexpr double rvv10_c_default = 0.0093;

// All density/length quantities use bohr; the combined energy convention is Ry.
// sigma = |grad n|^2. dq0_dsigma is NOT the coefficient 2*dq0/dsigma.
struct LocalValues
{
    double q0;
    double dq0_dn;
    double dq0_dsigma;
    double weight; // D = n/kappa_Ry^(3/2)
    double dweight_dn;
};

class LocalModel
{
  public:
    // Parameters are explicit: their choice belongs to the base-XC functional.
    LocalModel(double b, double c);
    LocalValues evaluate(double density, double sigma) const;
    double beta() const;

  private:
    double k_prefactor_;
    double gap_prefactor_;
    double weight_prefactor_;
    double beta_;
};

const std::array<double, channel_count>& q_mesh();

// The revised universal kernel paired with D_Ry at BOTH integration points.
double real_kernel(double q1, double q2, double distance);

struct Channel
{
    double theta;
    double dtheta_dn;
    double dtheta_dsigma;
};

struct KernelValues
{
    // Row-major symmetric matrices; derivative is with respect to |G| in bohr^-1.
    std::array<double, channel_count * channel_count> value;
    std::array<double, channel_count * channel_count> derivative;
};

class KernelTable
{
  public:
    KernelTable();
    KernelValues evaluate(double g) const;
    // SCF energy/potential need the kernel but not its |G| derivative.
    std::array<double, channel_count * channel_count> values(double g) const;
    double maximum_g() const;

  private:
    void interpolate(double g, double* value, double* derivative) const;
    static constexpr int pair_count = channel_count * (channel_count + 1) / 2;
    ModuleBase::CubicSpline spline_;
    std::array<int, pair_count> indices_;
};

class SplineBasis
{
  public:
    SplineBasis();
    Channel channel(const LocalValues& local, int index) const;

  private:
    // Owns the existing ABACUS spline; hence the complete type in this header.
    ModuleBase::CubicSpline spline_;
};
} // namespace Rvv10

#endif
