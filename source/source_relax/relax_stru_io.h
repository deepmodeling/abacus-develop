#ifndef RELAX_STRU_IO_H
#define RELAX_STRU_IO_H

#include "source_base/matrix.h"
#include "source_cell/unitcell.h"
#include "source_io/module_parameter/input_parameter.h"

#include <string>

/**
 * @brief Helpers for writing structure (STRU/CIF) files during relaxation.
 *
 * These free functions extract the file-writing logic shared by
 * Relax_Driver::stru_out and Relax_Driver::final_out so that the two
 * call sites only orchestrate when and which files to write.
 */
namespace relax_stru_io
{
    /**
     * @brief Build the header comment for a structure file.
     *
     * Contains the ABACUS version, a timestamp, the relaxation step number,
     * the total energy in eV and the 3x3 stress tensor in kbar.
     *
     * @param istep Current (zero-based) relaxation step; printed as istep + 1.
     * @param etot Total energy in Ry.
     * @param stress Stress tensor (3x3) in Ry/Bohr^3.
     * @param is_final If true, the step label is marked "(FINAL)".
     * @return The formatted header string.
     */
    std::string build_stru_header(const int istep,
                                  const double etot,
                                  const ModuleBase::matrix& stress,
                                  const bool is_final);

    /**
     * @brief Whether orbital output is required in the structure file.
     *
     * True for lcao / lcao_in_pw bases, and for pw with nao-initialized
     * wavefunctions.
     */
    bool need_orbital(const Input_para& inp);

    /**
     * @brief Write one structure file in STRU (out_stru == 1) or CIF
     *        (out_stru == 2) format.
     *
     * @param ucell Unit cell to write.
     * @param inp Input parameters (nspin, out_mul, calculation).
     * @param filename Output file path.
     * @param header Header comment built by build_stru_header.
     * @param force Atomic forces (used for STRU output).
     * @param need_orb Whether to write orbitals (from need_orbital).
     * @param deepks_setorb DeePKS setorb flag forwarded to the STRU writer.
     * @param my_rank MPI rank of the calling process.
     */
    void write_stru(UnitCell& ucell,
                    const Input_para& inp,
                    const std::string& filename,
                    const std::string& header,
                    const ModuleBase::matrix& force,
                    const bool need_orb,
                    const bool deepks_setorb,
                    const int my_rank);
} // namespace relax_stru_io

#endif
