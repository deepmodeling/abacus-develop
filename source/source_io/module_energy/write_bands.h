#ifndef WRITE_BANDS_H
#define WRITE_BANDS_H

#include "source_io/module_parameter/input_parameter.h"
#include "source_base/matrix.h"
#include "source_cell/klist.h"

#include <fstream>

namespace ModuleIO
{


/**
 * @brief write the eigenvalues of every spin channel to band.txt (nspin=1/4)
 *        or to bands1.txt/bands2.txt (nspin=2).
 *
 * @param inp input parameters, out_band[0] enables the output and out_band[1]
 *            sets the number of printed decimals
 * @param ekb eigenvalues of k points and bands
 * @param kv klist
 * @param ofs_running running log stream
 * @param istep zero-based index of the current ionic step, printed one-based
 *              in the per-step provenance header
 */
void write_bands(const Input_para& inp,
    const ModuleBase::matrix& ekb,
    const K_Vectors& kv,
    std::ofstream& ofs_running,
    const int istep);

/**
 * @brief calculate the band structure
 *
 * @param is spin index is = 0 or 1
 * @param eig_file file to save the band structure
 * @param nband number of bands
 * @param fermie fermi energy
 * @param precision precision of the output
 * @param ekb eigenvalues of k points and bands
 * @param kv klist
 * @param header provenance header written once in front of the k-point block;
 *               pass an empty string to omit it
 */
void nscf_bands(const int& is,
               const std::string &eig_file,
               const int& nband,
               const double& fermie,
               const int& precision,
               const ModuleBase::matrix& ekb,
               const K_Vectors& kv,
               const std::string& header);
}

#endif
