#ifndef ABACUS_SMOOTH_CUT_COULOMB_H
#define ABACUS_SMOOTH_CUT_COULOMB_H

#include <array>

namespace Conv_Coulomb_Pot_K
{
// k^2 v(k)/(4*pi) for f(r)/r, f=erfc(log(r/Rc)/log(width))/2.
double log_erfc_cut_factor(double wave_number, double cutoff_radius, double width_factor);

// Half the shortest nonzero superlattice vector; row vectors are in Bohr.
double bvk_inscribed_radius(const std::array<std::array<double, 3>, 3>& bvk_lattice);
}

#endif
