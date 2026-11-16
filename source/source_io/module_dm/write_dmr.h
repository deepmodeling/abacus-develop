#ifndef MODULE_IO_WRITE_DMR_H
#define MODULE_IO_WRITE_DMR_H

#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_hamilt/module_hcontainer/hcontainer.h"
#include "source_cell/unitcell.h"

#include <fstream>
#include <string>

namespace ModuleIO
{

/**
 * Generates a filename for the DMR output based on the given parameters.
 *
 * @param out_type The output type. 1: csr, 2: npz.
 * @param ispin    The spin value, starting from 0.
 * @param numbered Whether the ionic step is embedded in the filename as a g{step} suffix
 *                 (true for out_app_flag == 2).
 * @param istep    The ION step (default: -1), starting from 0.
 * @return         The generated filename as a string.
 */
std::string dmr_gen_fname(const int out_type,
    const int ispin,
    const bool numbered = false,
    const int istep = -1);

/**
 * Writes HContainer to a csr file.
 *
 * @param fname The name of the file to write the CSR representation to.
 * @param ucell stands for unitcell
 * @param precision Control the ouptut precision
 * @param dm_serial A pointer to the Hamiltonian container.
 * @param istep The current step number.
 * @param append Whether to append to the existing file. The first written step
 *               (istep <= 0) always truncates.
 */
void write_dmr_csr(std::string& fname,
        const UnitCell *ucell,
        const int precision,
        hamilt::HContainer<double>* dm_serial,
        const int istep,
        const int ispin,
        const int nspin,
        const bool append);

/**
 * Writes DMR to a file.
 *
 * @param dmr The 2D block parallel matrix representing the density matrix. The first dimension is the spin index.
 * @param ucell stands for unitcell
 * @param precision Control the output precision
 * @param paraV The parallel 2D object.
 * @param out_app_flag Output mode: 0 = fixed filename, overwritten at every write; 1 = fixed filename,
 *                    appended; 2 = per-step g{step} filename.
 * @param iat2iwt The mapping from atom index to the first local orbital index.
 * @param nat The number of atoms.
 * @param dmr_dir The directory the DM(R) files are written into.
 * @param istep The ION step, starting from 0.
 */
void write_dmr(const std::vector<hamilt::HContainer<double>*> dmr,
               const UnitCell *ucell,
               const int precision,
               const Parallel_Blkcyclic& paraV,
               const int out_app_flag,
               const int* iat2iwt,
               const int nat,
               const std::string& dmr_dir,
               const int istep,
               std::ofstream& ofs_running);
} // namespace ModuleIO

#endif
