#ifndef WRITE_EIG_OCC_H
#define WRITE_EIG_OCC_H
#include "source_base/matrix.h"
#include "source_cell/klist.h"

#include <string>

namespace ModuleIO
{
    /// @brief print eigenvalues and occupations to the running log
    /// @param nbands number of bands requested by INPUT
    /// @param nspin number of spin channels
    void write_eig_iter(const ModuleBase::matrix &ekb,
        const ModuleBase::matrix &wg,
        const K_Vectors& kv,
        const int nbands,
        const int nspin);

    /// @brief write eigenvalues and occupations to <out_dir>/eig_occ.txt
    /// @param nbands number of bands requested by INPUT
    /// @param nspin number of spin channels
    /// @param out_dir directory the file is written into
    void write_eig_file(const ModuleBase::matrix &ekb,
            const ModuleBase::matrix &wg,
            const K_Vectors& kv,
            const int nbands,
            const int nspin,
            const std::string& out_dir,
            const int istep);
}

#endif
