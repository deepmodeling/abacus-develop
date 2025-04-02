#ifndef SPARSE_FORMAT_EXX_H
#define SPARSE_FORMAT_EXX_H

#ifdef __EXX

#include <RI/global/Global_Func-2.h>
#include <RI/ri/Cell_Nearest.h>
#include <array>
#include <map>
#include <vector>

#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_HS_arrays.hpp"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_cell/unitcell.h"
namespace sparse_format
{

template <typename Tdata>
void cal_HR_exx(const UnitCell& ucell,
    const Parallel_Orbitals& pv,
    LCAO_HS_Arrays& HS_Arrays,
    const int& current_spin,
    const double& sparse_thr,
    const int (&nmp)[3],
    const std::vector<std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<Tdata>>>>& Hexxs);

}
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_hamilt.hpp"
#endif
#endif
