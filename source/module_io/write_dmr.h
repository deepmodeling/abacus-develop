#ifndef MODULE_IO_WRITE_DMR_H
#define MODULE_IO_WRITE_DMR_H
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"

#include <string>

namespace ModuleIO
{

/**
 * Generates a filename for the DMR output based on the given parameters.
 *
 * @param out_type The output type. 1: csr, 2: npz.
 * @param ispin    The spin value, starting from 0.
 * @param append   A boolean indicating whether append the data to one file or create a new file.
 * @param istep    The ION step (default: -1), starting from 0.
 * @return         The generated filename as a string.
 */
std::string dmr_gen_fname(const int out_type, const int ispin, const bool append = true, const int istep = -1);

/**
 * Writes HContainer to a csr file.
 *
 * @param fname The name of the file to write the CSR representation to.
 * @param dm_serial A pointer to the Hamiltonian container.
 * @param istep The current step number.
 */
void write_dmr_csr(std::string& fname, hamilt::HContainer<double>* dm_serial, const int istep);

/**
 * Writes DMR to a file.
 *
 * @param dmr The 2D block parallel matrix representing the density matrix. The first dimension is the spin index.
 * @param paraV The parallel 2D object.
 * @param out_type The output file type. 1: csr, 2: npz.
 * @param sparse Whether output the sparse DM.
 * @param ispin The spin index, starting from 0.
 * @param append Whether to append the data to an existing file or create a new file. The file name is related to this
 * flag.
 * @param istep The ION step, starting from 0.
 */
void write_dmr(const std::vector<hamilt::HContainer<double>*> dmr,
               const Parallel_2D& paraV,
               const bool out_csr,
               const bool out_npz,
               const bool append,
               const int* iat2iwt,
               const int* nat,
               const int istep);
} // namespace ModuleIO

#endif