#ifndef XC_NCGGA_RADIAL_H
#define XC_NCGGA_RADIAL_H

#include <array>

namespace ModuleXC
{

struct NcggaRadialPoint
{
    // Value, gradient, and Hessian of the same radial scalar map.  At zero,
    // direction is represented by the zero vector because the scalar map has
    // a unique zero gradient and zero Hessian there.
    double value = 0.0;
    std::array<double, 3> direction = {{0.0, 0.0, 0.0}};
    std::array<double, 3> gradient = {{0.0, 0.0, 0.0}};
    double transverse_hessian = 0.0;
    double radial_hessian = 0.0;

    double jacobian(const int row, const int column) const;
};

// For r = |magnetization| and x = r / eta, the returned scalar is
//   eta * x^3 * (3 x^2 - 8 x + 6), r < eta,
//   r,                                  r >= eta.
// The splice is C2 at both r = 0 and r = eta.  Eta is explicit so this
// mathematical primitive does not choose policy for any XC mode.
NcggaRadialPoint make_ncgga_radial_point(const std::array<double, 3>& magnetization, const double eta);

// The C2 regularization scale is part of the gga_grad=2 LCA functional, not a
// divide-by-zero guard.  Other noncollinear modes do not inherit this policy.
double ncgga_lca_radial_eta();

struct NcggaSpinMapPoint
{
    NcggaRadialPoint radial;
    double absolute_density = 0.0;
    double clipped_magnitude = 0.0;
    std::array<double, 2> spin_density = {{0.0, 0.0}};
    double density_sign = 0.0;
    bool saturated = true;

    // spin is 0/1 for up/down; channel is 0 for the raw total density and
    // 1..3 for mx,my,mz.
    double jacobian(const int spin, const int channel) const;
};

// Compose a raw total density t=n+rho_core with one radial magnetization map:
//   a       = |t|,
//   c       = min(radial.value, a),
//   rho_up  = (a+c)/2,
//   rho_down= (a-c)/2.
// The Jacobian follows that exact graph. At the abs kink t=0 it selects zero;
// at the clipping kink radial.value=a it selects the saturated branch. Eta and
// the radial policy remain explicit in make_ncgga_radial_point.
NcggaSpinMapPoint make_ncgga_spin_map_point(const double total_density, const NcggaRadialPoint& radial);

} // namespace ModuleXC

#endif
