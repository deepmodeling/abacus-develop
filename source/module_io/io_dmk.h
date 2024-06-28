#ifndef DM_IO_H
#define DM_IO_H

#include "module_basis/module_ao/parallel_2d.h"
#include "module_cell/unitcell.h"

#include <string>

namespace ModuleIO
{
/**
 * @brief Reads the DMK file and stores the data in the provided arrays.
 *
 * @param nnrg The number of energy points.
 * @param trace_lo An array containing the lower trace indices.
 * @param gamma_only_local A boolean indicating whether to read only the gamma-only part of the DMK file.
 * @param nlocal The number of local orbitals.
 * @param nspin The number of spin components.
 * @param is The index of the spin component.
 * @param fn The filename of the DMK file.
 * @param DM A 3D array to store the DMK data.
 * @param DM_R A 2D array to store the real part of the DMK data.
 * @param ef The Fermi energy.
 * @param ucell A pointer to the UnitCell object.
 */
void read_dmk(
#ifdef __MPI
    const int nnrg,
    const int* trace_lo,
#endif
    const bool gamma_only_local,
    const int nlocal,
    const int nspin,
    const int& is,
    const std::string& fn,
    double*** DM,
    double** DM_R,
    double& ef,
    const UnitCell* ucell);

/**
 * @brief Generates the filename for the DMK file based on the given parameters.
 *
 * @param gamma_only A boolean indicating whether to generate the filename for the gamma-only part.
 * @param ispin The index of the spin component.
 * @param ik The index of the k-point.
 * @return The generated filename.
 */
std::string dmk_gen_fname(const bool gamma_only, const int ispin, const int ik);

/**
 * @brief Writes the DMK data to a file.
 *
 * @tparam T The type of the DMK data.
 * @param dmk A vector containing the DMK data. The first dimension is nspin*nk, and the second dimension is
 * nlocal*nlocal. DMK is parallel in 2d-block type if using MPI.
 * @param precision The precision of the output of DMK.
 * @param efs A vector containing the Fermi energies, and should have the same size as the number of SPIN.
 * @param ucell A pointer to the UnitCell object.
 * @param pv The Parallel_2D object. The 2d-block parallel information of DMK.
 */
template <typename T>
void write_dmk(const std::vector<std::vector<T>> dmk,
               const int precision,
               const std::vector<double> efs,
               const UnitCell* ucell,
               const Parallel_2D& pv);
} // namespace ModuleIO
#endif
