#ifndef MODULE_IO_WRITE_DMR_H
#define MODULE_IO_WRITE_DMR_H
#include <string>

#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_basis/module_ao/parallel_orbitals.h"

namespace ModuleIO
{

/**
 * Generates a filename for the DMR output based on the given parameters.
 *
 * @param out_type The output type. 1: csr, 2: npz.
 * @param sparse   A boolean indicating whether the output is sparse.
 * @param ispin    The spin value, starting from 0.
 * @param append   A boolean indicating whether append the data to one file or create a new file.
 * @param istep    The ION step (default: -1), starting from 0.
 * @return         The generated filename as a string.
 */
std::string dmr_gen_fname(const int out_type, const bool sparse, const int ispin, const bool append=true, const int istep=-1);

/**
 * Writes DMR to a file in CSR format
 *
 * @param fname The output stream to write the CSR data to.
 * @param dm The Hamiltonian container to write.
 * @param istep The ION step, starting from 0.
 */
void write_dmr_csr(std::string& fname, const hamilt::HContainer<double>& dm, const int istep);

/**
 * Writes DMR to a file.
 *
 * @param dm The 2D block parallel matrix representing the density matrix.
 * @param out_type The output file type. 1: csr, 2: npz.
 * @param sparse Whether output the sparse DM.
 * @param ispin The spin index, starting from 0.
 * @param append Whether to append the data to an existing file or create a new file. The file name is related to this flag.
 * @param istep The ION step, starting from 0.
 * @param pv The Parallel_Orbitals object.
 */
void write_dmr(
    const hamilt::HContainer<double>& dm,
    const int out_type,
    const bool sparse,
    const int ispin,
    const bool append,
    const int istep, // start from 0
    const Parallel_Orbitals& pv);
}

#endif