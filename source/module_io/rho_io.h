#ifndef RHO_IO_H
#define RHO_IO_H
#include <string>
#include "module_cell/unitcell.h"
#ifdef __MPI
#include "module_hamilt_pw/hamilt_pwdft/parallel_grid.h"
#endif

namespace ModuleIO
{
bool read_rho(
#ifdef __MPI
    Parallel_Grid* Pgrid,
#endif
    int my_rank,
    std::string esolver_type,
    int rank_in_stogroup,
    const int& is,
    std::ofstream& ofs_running,
    const int& nspin,
    const std::string& fn,
    double* rho,
    int& nx,
    int& ny,
    int& nz,
    double& ef,
    const UnitCell* ucell,
    int& prenspin);
}

#endif
