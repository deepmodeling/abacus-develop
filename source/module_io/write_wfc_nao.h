#ifndef WRITE_WFC_NAO_H
#define WRITE_WFC_NAO_H
#include "module_base/matrix.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_psi/psi.h"
#include "module_base/vector3.h"
#include <vector>

namespace ModuleIO
{

/**
 * Generates the filename for the output of LCAO WFC.
 *
 * @param out_type The type of output.
 * @param gamma_only Whether gamma_only job.
 * @param out_app_flag Flag indicating whether to append to existing file.
 * @param ik The index of the k-point, and starting from 0.
 * @param istep The index of the ION step, and starting from 0. If < 0, the step number is not included in the file name.
 * @return The generated filename.
 */
std::string wfc_nao_gen_fname(const int out_type,
                               const bool gamma_only,
                               const bool out_app_flag,
                               const int ik,
                               const int istep=-1);

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
                    const psi::Psi<T>& psi,
                    const ModuleBase::matrix& ekb,
                    const ModuleBase::matrix& wg,
                    const std::vector<ModuleBase::Vector3<double>>& kvec_c,
                    const Parallel_Orbitals& pv,
                    const int istep=-1) ;

void wfc_nao_write2file(const std::string &name, const double* ctot, const int nlocal, const int ik, const ModuleBase::matrix& ekb, const ModuleBase::matrix& wg, bool writeBinary);
void wfc_nao_write2file_complex(const std::string &name, const std::complex<double>* ctot, const int nlocal,const int &ik, const ModuleBase::Vector3<double> &kvec_c, const ModuleBase::matrix& ekb, const ModuleBase::matrix& wg, bool writeBinary=false);
}// namespace ModuleIO
#endif