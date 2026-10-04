#ifndef VXC_OP_TOOLS_H
#define VXC_OP_TOOLS_H

#include "source_base/parallel_2d.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_base/module_container/base/third_party/blas.h"
#include "source_base/parallel_reduce.h"
#include "source_cell/klist.h"

#include <complex>
#include <string>
#include <vector>

namespace ModuleIO
{

/// @brief Set up a 2D parallel distribution for MO-space matrices (nbands x nbands)
/// @param pv Parallel_Orbitals describing the AO distribution
/// @param nbands Number of bands
/// @param p2d Output Parallel_2D object
void set_para2d_MO(const Parallel_Orbitals& pv, const int nbands, Parallel_2D& p2d);

/// @brief Extract real part (for complex numbers)
inline double get_real(const std::complex<double>& c)
{
    return c.real();
}

/// @brief Extract real part (for real numbers, identity)
inline double get_real(const double& d)
{
    return d;
}

/// @brief Compute c^dagger * V * c in MO basis
/// @tparam T Data type (double or std::complex<double>)
/// @param V AO-space matrix V (nbasis x nbasis)
/// @param c AO-to-MO transformation coefficients (nbasis x nbands)
/// @param nbasis Number of AO basis functions
/// @param nbands Number of bands
/// @param pv Parallel_Orbitals for AO distribution
/// @param p2d Parallel_2D for MO distribution
/// @return MO-space matrix c^dagger * V * c (nbands x nbands)
template <typename T>
std::vector<T> cVc(const T* V,
                   const T* c,
                   const int nbasis,
                   const int nbands,
                   const Parallel_Orbitals& pv,
                   const Parallel_2D& p2d);

/// @brief Extract diagonal elements (orbital energies) from MO matrix
/// @tparam T Data type (double or std::complex<double>)
/// @param ik K-point index (for occupancy weights)
/// @param nbands Number of bands
/// @param mat_mo MO-space matrix
/// @param p2d Parallel_2D distribution
/// @return Vector of orbital energies (length nbands)
template <typename T>
std::vector<double> orbital_energy(const int ik,
                                   const int nbands,
                                   const std::vector<T>& mat_mo,
                                   const Parallel_2D& p2d);

/// @brief Sum over bands with occupation weights
/// @tparam T Data type (double or std::complex<double>)
/// @param ik K-point index
/// @param mat_mo MO-space matrix
/// @param p2d Parallel_2D distribution
/// @param wg Occupation weights
/// @return Weighted sum of diagonal elements
template <typename T>
double all_band_energy(const int ik,
                       const std::vector<T>& mat_mo,
                       const Parallel_2D& p2d,
                       const ModuleBase::matrix& wg);

/// @brief Write orbital energies to file in LibRPA format
/// @param kv K-point list
/// @param nspin0 Number of independent spin channels (1 or 2)
/// @param nbands Number of bands
/// @param e_orb Orbital energies per k-point
/// @param term Term label (e.g., "vxc", "kinetic")
/// @param label Additional label (e.g., "local", "exx")
/// @param global_out_dir Output directory
/// @param app Append mode flag
void write_orb_energy(const K_Vectors& kv,
                      const int nspin0,
                      const int nbands,
                      const std::vector<std::vector<double>>& e_orb,
                      const std::string& term,
                      const std::string& label,
                      const std::string& global_out_dir,
                      const bool app = false);

} // namespace ModuleIO

#endif // VXC_OP_TOOLS_H
