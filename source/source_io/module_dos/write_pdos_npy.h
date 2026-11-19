#ifndef WRITE_PDOS_NPY_H
#define WRITE_PDOS_NPY_H

#include <string>
#include <vector>

class UnitCell;
namespace ModuleBase
{
class matrix;
}

namespace ModuleIO
{
std::string write_pdos_npy(const UnitCell& cell, const std::vector<ModuleBase::matrix>& pdos,
                           int nspin, double emin, double de, const std::string& out_dir,
                           int istep, int out_app_flag);
}
#endif
