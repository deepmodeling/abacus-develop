#ifndef MODULE_IO_TD_EIJ_IO_H
#define MODULE_IO_TD_EIJ_IO_H

#include <complex>
#include <ostream>

namespace ModuleIO
{
/** @brief Eij output settings and one-based identifiers for one electronic state. */
struct EijOutputInfo
{
    double threshold;
    int electronic_step;
    int scf_iteration;
    int kpoint;
    int nspin;
    int spin;
};

/** @brief Print selected elements of a column-major Hamiltonian in the orbital basis, in Ry. */
void print_eij(const std::complex<double>* eij,
               const int nband,
               const int leading_dimension,
               const EijOutputInfo& output_info,
               std::ostream& ofs);
} // namespace ModuleIO
#endif
