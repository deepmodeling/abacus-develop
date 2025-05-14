#ifndef DM_IO_H
#define DM_IO_H

#include "module_base/parallel_2d.h"
#include "module_cell/unitcell.h"

#include <string>
#include <vector>

namespace ModuleIO {

/**
 * @brief Generates the filename for the DMK file based on the given parameters.
 *
 * @param gamma_only Whether the calculation is gamma_only.
 * @param ispin The index of the spin component.
 * @param ik The index of the k-point.
 * @return The generated filename.
 */
std::string dmk_gen_fname(const bool gamma_only, const int ispin, const int ik);

/**
 * @brief Writes the unit cell information to a DMK file.
 *
 * @param ofs The output file stream.
 * @param ucell A pointer to the UnitCell object.
 */
void dmk_write_ucell(std::ofstream& ofs, const UnitCell* ucell);

/**
 * @brief Reads the unit cell information lines in a DMK file.
 *
 * @param ifs The input file stream.
 */
void dmk_read_ucell(std::ifstream& ifs);

/**
 * @brief Read one double from a file.
 */
void dmk_readData(std::ifstream& ifs, double& data);

/**
 * @brief Read one complex double from a file.
 */
void dmk_readData(std::ifstream& ifs, std::complex<double>& data);

/**
 * @brief Reads the DMK data from a file.
 *
 * @tparam T The type of the DMK data.
 * @param nspin The number of spin components.
 * @param nk The number of k-points.
 * @param pv The Parallel_2D object. Will get the global size and local size
 * from it, and seperate the data into different processors accordingly.
 * @param dmk_dir The directory path of the DMK file.
 * @param dmk A vector to store the DMK data. If use MPI parallel, the data will
 * be seperated into different processors based on the Parallel_2D object.
 * @return True if the DMK data is successfully read, false otherwise.
 */
template <typename T>
bool read_dmk(const int nspin,
              const int nk,
              const Parallel_2D& pv,
              const std::string& dmk_dir,
              std::vector<std::vector<T>>& dmk);

/**
 * @brief Writes the DMK data to a file.
 *
 * @tparam T The type of the DMK data.
 * @param dmk A vector containing the DMK data. The first dimension is nspin*nk,
 * and the second dimension is nlocal*nlocal. DMK is parallel in 2d-block type
 * if using MPI.
 * @param precision The precision of the output of DMK.
 * @param efs A vector containing the Fermi energies, and should have the same
 * size as the number of SPIN.
 * @param ucell A pointer to the UnitCell object.
 * @param pv The Parallel_2D object. The 2d-block parallel information of DMK.
 */
template <typename T>
void write_dmk(const std::vector<std::vector<T>>& dmk,
               const int precision,
               const std::vector<double>& efs,
               const UnitCell* ucell,
               const Parallel_2D& pv);

} // namespace ModuleIO

#endif // IO_DMK_H
