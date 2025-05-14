#ifndef READ_WFC_LCAO_H
#define READ_WFC_LCAO_H

// serial
#include "module_base/vector3.h"

#include <complex>
#include <string>
#include <vector>
#ifdef __MPI
// parallelization
#include "module_base/scalapack_connector.h"
#include "module_base/parallel_2d.h"
#endif

/**
 * @brief This class has two functions: restart psi from the previous calculation, and write psi to the disk.
 *
 */
namespace ModuleIO
{
// only when you know why you need the T, you can write function with template parameter,
// otherwise, you should overload the function for different types
// For example in this case, ONLY lowf support to be std::complex<double> and std::complex<float>,
// not ekb, occ, wk and kvec_c.
/**
 * @brief Read the wavefunction coefficients from the file (for complex wavefunction coefficients)
 *
 * @tparam T
 * @param flowf [in] file name like "LOWF_K_*.txt", dumped from ABACUS INPUT out_wfc_lcao 1
 * @param ik [out] the index of k points
 * @param kvec_c [out] the k vector in Cartesian coordinates
 * @param nbands [out] the number of bands
 * @param nbasis [out] the number of orbitals
 * @param lowf [out] wavefunction coefficients
 * @param ekb [out] eigenvalues
 * @param occ [out] occupations
 * @param wk [out] weight of k points
 */
template <typename T>
void read_abacus_lowf(const std::string& flowf, int& ik, ModuleBase::Vector3<double>& kvec_c, int& nbands, int& nbasis,
                      std::vector<std::complex<T>>& lowf, std::vector<double>& ekb, std::vector<double>& occ,
                      double& wk);
/**
 * @brief Read the wavefunction coefficients from the file (for real wavefunction coefficients)
 *
 * @tparam T
 * @param flowf [in] file name like "LOWF_K_*.txt", dumped from ABACUS INPUT out_wfc_lcao 1
 * @param ik [out] the index of k points
 * @param kvec_c [out] the k vector in Cartesian coordinates
 * @param nbands [out] the number of bands
 * @param nbasis [out] the number of orbitals
 * @param lowf [out] wavefunction coefficients
 * @param ekb [out] eigenvalues
 * @param occ [out] occupations
 * @param wk [out] weight of k points
 */
template <typename T>
void read_abacus_lowf(const std::string& flowf, int& ik, ModuleBase::Vector3<double>& kvec_c, int& nbands, int& nbasis,
                      std::vector<T>& lowf, std::vector<double>& ekb, std::vector<double>& occ, double& wk);
// the two functions above will return nbands, nbasis, lowf, ekb, occ and wk.
// the lowf is actually lowf_glb, which means the global matrix (ScaLAPACK convention), need to distribute
// to the local matrix (2D-block-cyclic parallel distribution) in the following function.

// only-MPI-visible function, because the use of comm_world
#ifdef __MPI
/**
 * @brief Restart the wavefunction coefficients from the file (MPI 2D-BCD version)
 *
 * @tparam T: datatype of the wavefunction coefficients, can be double, float, std::complex<double> or
 * std::complex<float>
 * @param out_dir [in] the directory where the wavefunction coefficients are stored
 * @param p2d [in] the 2D parallel distribution
 * @param nks [in] the number of k points
 * @param nbands [out] the number of bands
 * @param nbasis [out] the number of orbitals
 * @param lowf_loc [out] the local wavefunction coefficients, can be used to construct psi, see constructor No.8
 * @param ekb [out] the eigenvalues
 * @param occ [out] the occupations
 * @param kvec_c [out] the k vectors in Cartesian coordinates
 * @param wk [out] the weight of k points
 *
 * @warning Cpxgemr2d not implemented yet
 */
template <typename T>
void restart_from_file(const std::string& out_dir, // hard-code the file name to be LOWF_K_*.txt?
                       const Parallel_2D& p2d, const int& nks, int& nbands, int& nbasis, std::vector<T>& lowf_loc,
                       std::vector<double>& ekb, std::vector<double>& occ,
                       std::vector<ModuleBase::Vector3<double>>& kvec_c, std::vector<double>& wk);
#endif
// serial version, can always present
/**
 * @brief Restart the wavefunction coefficients from the file (serial version)
 *
 * @tparam T: datatype of the wavefunction coefficients, can be double, float, std::complex<double> or
 * std::complex<float>
 * @param out_dir [in] the directory where the wavefunction coefficients are stored
 * @param nks [in] the number of k points
 * @param nbands [out] the number of bands
 * @param nbasis [out] the number of orbitals
 * @param lowf_loc [out] the local wavefunction coefficients, can be used to construct psi, see constructor No.8
 * @param ekb [out] the eigenvalues
 * @param occ [out] the occupations
 * @param kvec_c [out] the k vectors in Cartesian coordinates
 * @param wk [out] the weight of k points
 */
template <typename T>
void restart_from_file(const std::string& out_dir, // hard-code the file name to be LOWF_K_*.txt?
                       const int& nks, int& nbands, int& nbasis, std::vector<T>& lowf, std::vector<double>& ekb,
                       std::vector<double>& occ, std::vector<ModuleBase::Vector3<double>>& kvec_c,
                       std::vector<double>& wk);
} // namespace ModuleIO
#endif