#ifndef READ_WFC_NAO_H
#define READ_WFC_NAO_H

#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_psi/psi.h"
#include "source_base/matrix.h"

// mohan add 2010-09-09
namespace ModuleIO
{
/**
 * @brief Reads a single data value from an input file stream.
 * 
 * @param ifs The input file stream to read from.
 * @param data The variable to store the read data value.
 */
void read_wfc_nao_one_data(std::ifstream& ifs, float& data);
void read_wfc_nao_one_data(std::ifstream& ifs, double& data);

/**
 * @brief Reads a single complex data value from an input file stream.
 * 
 * @param ifs The input file stream to read from.
 * @param data The variable to store the read complex data value.
 */
void read_wfc_nao_one_data(std::ifstream& ifs, std::complex<float>& data);
void read_wfc_nao_one_data(std::ifstream& ifs, std::complex<double>& data);

/**
 * @brief Reads the wavefunction coefficients from an input file.
 * 
 * @tparam T The type of the wavefunction coefficients.
 * @param global_readin_dir The global directory for reading input files.
 * @param ParaV The parallel orbitals object.
 * @param psid The Psi object to store the wavefunction coefficients.
 * @param pelec Pointer to the ElecState object.
 * @param skip_band From which band to start reading.
 * @return True if the wavefunction coefficients are successfully read, false otherwise.
 */
template <typename T>
bool read_wfc_nao(
    const std::string& global_readin_dir,
    const Parallel_Orbitals& ParaV,
    psi::Psi<T>& psid,
	ModuleBase::matrix& ekb,
    ModuleBase::matrix& wg,
	const std::vector<int> &ik2iktot,
	const int nkstot,
	const int nspin,
    const int skip_band = 0,
    const int istep = -1);

} // namespace ModuleIO

#endif
