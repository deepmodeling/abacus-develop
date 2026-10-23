#ifndef WRITE_WFC_NAO_H
#define WRITE_WFC_NAO_H
#include "source_base/matrix.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_psi/psi.h"
#include "source_base/vector3.h"
#include <string>
#include <vector>

namespace ModuleIO
{

    /**
 * Writes the wavefunction coefficients for the LCAO method to a file.
 * Will loop all k-points by psi.get_nk().
 * The nbands are determined by ekb.nc.
 * The nlocal is determined by psi.get_nbasis() if not compiled with MPI, otherwise it is determined by pv->desc[2].
 *
 * @param out_type The output file type. 1 for text, 2 for binary.
 * @param psi The Psi object containing the wavefunction coefficients.
 * @param ekb The matrix of Kohn-Sham eigenvalues.
 * @param wg The matrix of Kohn-Sham eigenvectors.
 * @param kvec_c The vector of k-points in Cartesian coordinates.
 * @param pv The Parallel_Orbitals object containing additional information.
 * @param istep The current step number. if < 0, the step number is not included in the file name.
 */
template <typename T>
void write_wfc_nao(const int out_type,
            const int out_app_flag,
            const psi::Psi<T>& psi,
            const ModuleBase::matrix& ekb,
            const ModuleBase::matrix& wg,
            const std::vector<ModuleBase::Vector3<double>>& kvec_c,
            const std::vector<int> &ik2iktot,
            const int nkstot,
            const Parallel_Orbitals& pv,
            const int nspin,
            const int istep=-1) ;

/**
 * @brief write real-valued Gamma-only wavefunction coefficients to a file.
 *
 * @param name output file name
 * @param ctot global wavefunction coefficient array, size nbands*nlocal
 * @param nlocal number of orbitals (basis functions)
 * @param ik local k-point index, used to index ekb/wg rows
 * @param ekb eigenvalue matrix
 * @param wg occupation matrix
 * @param writeBinary true for the binary format, false for text
 * @param append_flag true to append to an existing file, false to truncate
 * @param header provenance comment block written before the data in text
 *        mode; ignored in binary mode
 */
void wfc_nao_write2file(const std::string& name,
                        const double* ctot,
                        const int nlocal,
                        const int ik,
                        const ModuleBase::matrix& ekb,
                        const ModuleBase::matrix& wg,
                        const bool& writeBinary,
                        const bool& append_flag,
                        const std::string& header);

/**
 * @brief write complex wavefunction coefficients to a file.
 *
 * @param name output file name
 * @param ctot global wavefunction coefficient array, size nbands*nlocal
 * @param nlocal number of orbitals (basis functions)
 * @param ik local k-point index, used to index ekb/wg rows
 * @param kvec_c Cartesian k-vector of the k-point
 * @param ekb eigenvalue matrix
 * @param wg occupation matrix
 * @param writeBinary true for the binary format, false for text
 * @param append_flag true to append to an existing file, false to truncate
 * @param header provenance comment block written before the data in text
 *        mode; ignored in binary mode
 */
void wfc_nao_write2file_complex(const std::string& name,
                                const std::complex<double>* ctot,
                                const int nlocal,
                                const int& ik,
                                const ModuleBase::Vector3<double>& kvec_c,
                                const ModuleBase::matrix& ekb,
                                const ModuleBase::matrix& wg,
                                const bool& writeBinary,
                                const bool& append_flag,
                                const std::string& header);
}// namespace ModuleIO
#endif
