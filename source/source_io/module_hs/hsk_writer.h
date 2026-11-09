#ifndef HSK_WRITER_H
#define HSK_WRITER_H

#include <fstream>
#include <string>
#include <vector>

#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_cell/klist.h"
#include "source_cell/unitcell.h"
#include "source_hamilt/hamilt.h"

// mohan note: this file holds the H(k)/S(k) writer
// (formerly write_hsk in write_hs.h / write_hs.hpp).

namespace ModuleIO
{

/**
 * @brief Write H(k) and S(k) matrices for each k-point to disk.
 *
 * Text output format is aligned with write_dmk: full metadata block
 * (nspin, ik, k-point coords, weight, Fermi energy, basis count),
 * followed by a UnitCell lattice block and then a flat matrix body
 * (4 complex or 8 double values per line). Binary output continues
 * to delegate to save_mat so its layout stays unchanged for callers
 * that only consume binary files.
 *
 * @param kv          K-point container (coords, weights, counts).
 * @param ucell       Pointer to the UnitCell for the lattice block;
 *                    may be nullptr if the caller cannot provide one.
 * @param efs         Fermi energies — one per spin channel.
 */
template <typename T>
void write_hsk(const std::string& global_out_dir,
               const int nspin,
               const int nks,
               const int nkstot,
               const std::vector<int>& ik2iktot,
               const std::vector<int>& isk,
               hamilt::Hamilt<T>* p_hamilt,
               const Parallel_Orbitals& pv,
               const K_Vectors& kv,
               const UnitCell* ucell,
               const std::vector<double>& efs,
               const bool gamma_only,
               const int out_app_flag,
               const int istep,
               const int out_type,
               const int precision,
               const int nlocal,
               const std::string& ks_solver,
               const int drank,
               std::ofstream& ofs_running);

} // namespace ModuleIO

#endif // HSK_WRITER_H
