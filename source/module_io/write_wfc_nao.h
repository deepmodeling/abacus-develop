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
			const bool out_app_flag,
			const psi::Psi<T>& psi,
			const ModuleBase::matrix& ekb,
			const ModuleBase::matrix& wg,
			const std::vector<ModuleBase::Vector3<double>>& kvec_c,
			const std::vector<int> &ik2iktot,
			const int nkstot,
			const Parallel_Orbitals& pv,
			const int nspin,
			const int istep=-1) ;

void wfc_nao_write2file(const std::string& name,
                        const double* ctot,
                        const int nlocal,
                        const int ik,
                        const ModuleBase::matrix& ekb,
                        const ModuleBase::matrix& wg,
                        const bool& writeBinary,
                        const bool& append_flag = false);

void wfc_nao_write2file_complex(const std::string& name,
                                const std::complex<double>* ctot,
                                const int nlocal,
                                const int& ik,
                                const ModuleBase::Vector3<double>& kvec_c,
                                const ModuleBase::matrix& ekb,
                                const ModuleBase::matrix& wg,
                                const bool& writeBinary = false,
                                const bool& append_flag = false);
}// namespace ModuleIO
#endif
