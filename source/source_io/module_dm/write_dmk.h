#ifndef WRITE_DMK_H
#define WRITE_DMK_H

#include "source_base/parallel_blkcyc.h"
#include "source_cell/unitcell.h"
#include "source_cell/klist.h"

#include <fstream>
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
std::string dmk_gen_fname(const bool gamma_only, const int ispin, const int nspin, const int ik, const int istep);

/**
 * @brief Read one double from a file.
 * @return false when the extraction fails (truncated or malformed body).
 */
bool dmk_readData(std::ifstream& ifs, double& data);

/**
 * @brief Read one complex double from a file.
 * @return false when the extraction fails (truncated, empty token or a token
 * that is not of the "(re,im)" form).
 */
bool dmk_readData(std::ifstream& ifs, std::complex<double>& data);

/**
 * @brief Reads the DMK data from a file.
 *
 * @tparam T The type of the DMK data.
 * @param nspin The number of spin components.
 * @param nk The number of k-points.
 * @param k-vector information.
 * @param pv The Parallel_Blkcyclic object. Will get the global size and local size
 * from it, and seperate the data into different processors accordingly.
 * @param column_major Whether the local DMK blocks must be stored column-major
 * (genelpa/elpa/scalapack_gvx and the other column-major solvers). The file body
 * is row-major, so this selects the destination index of each element.
 * @param dmk_dir The directory path of the DMK file.
 * @param istep The ionic step to read. With out_app_flag == 1 it selects the
 * frame of the appended multi-step file, and a negative value keeps the first
 * frame; with out_app_flag == 2 it has to match the g{istep+1} file name and
 * the frame marker it contains.
 * @param out_app_flag Must mirror the value used by write_dmk, because it
 * decides whether the step is part of the filename.
 * @param dmk A vector to store the DMK data. If use MPI parallel, the data will
 * be seperated into different processors based on the Parallel_Blkcyclic object.
 * @return True if the DMK data is successfully read, false otherwise.
 */
template <typename T>
bool read_dmk(const int nspin,
              const int nk,
              const K_Vectors &kv,
              const Parallel_Blkcyclic& pv,
              const bool column_major,
              const std::string& dmk_dir,
              const int istep,
              const int out_app_flag,
              std::vector<std::vector<T>>& dmk,
              std::ofstream &ofs_running);

/**
 * @brief Writes the DMK data to a file.
 *
 * @details The written matrix element is
 *          D(i,j) = sum_n w_n * conj(c_n,i) * c_n,j
 *          (complex conjugate on the first band index), matching the
 *          convention built by dm_from_psi; this is the complex conjugate of
 *          the physical one-body density matrix sum_n w_n c_n,i conj(c_n,j).
 *          The same line is also emitted as a '#' comment in the file header.
 *
 * @tparam T The type of the DMK data.
 * @param dmk A vector containing the DMK data. The first dimension is nspin*nk,
 * and the second dimension is nlocal*nlocal. DMK is parallel in 2d-block type
 * if using MPI.
 * @param k-vector information.
 * @param precision The precision of the output of DMK.
 * @param efs A vector containing the Fermi energies, and should have the same
 * size as the number of SPIN.
 * @param ucell A pointer to the UnitCell object.
 * @param pv The Parallel_Blkcyclic object. The 2d-block parallel information of DMK.
 * @param column_major Whether the local DMK blocks are stored column-major
 * (genelpa/elpa/scalapack_gvx and the other column-major solvers). The file
 * body is always row-major, so this selects the traversal order.
 * @param dmk_dir The directory the DMK files are written into.
 * @param out_app_flag Append-mode selector: 1 reuses the fixed filename
 * across ionic steps and appends a new record per step in multi-step
 * runs (MD or relax), 0 or 2 inserts the geometry-step index g{N}.
 */
template <typename T>
void write_dmk(const std::vector<std::vector<T>>& dmk,
               const K_Vectors &kv,
               const int precision,
               const std::vector<double>& efs,
               const UnitCell* ucell,
               const Parallel_Blkcyclic& pv,
               const bool column_major,
               const std::string& dmk_dir,
               const int istep,
               const int out_app_flag,
               std::ofstream& ofs_running);

} // namespace ModuleIO

#endif
