#ifndef NSCF_BAND_H
#define NSCF_BAND_H
#include "module_base/matrix.h"
#include "module_cell/klist.h"
#include "module_cell/parallel_kpoints.h"

namespace ModuleIO
{
/**
 * @brief calculate the band structure
 *
 * @param is spin index is = 0 or 1
 * @param out_band_dir directory to save the band structure
 * @param nband number of bands
 * @param fermie fermi energy
 * @param precision precision of the output
 * @param ekb eigenvalues of k points and bands
 * @param kv klist
 * @param Pkpoints parallel kpoints
 */
void nscf_band(const int& is,
               const std::string& out_band_dir,
               const int& nband,
               const double& fermie,
               const int& precision,
               const ModuleBase::matrix& ekb,
               const K_Vectors& kv,
               const Parallel_Kpoints* Pkpoints);
}

#endif
