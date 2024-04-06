#ifndef SPARSE_FORMAT_EXX_H 
#define SPARSE_FORMAT_EXX_H

#ifdef __EXX

#include <vector>
#include <map>
#include <array>

#include <RI/global/Global_Func-2.h>
#include <RI/ri/Cell_Nearest.h>

// use LCAO_Matrix
#include "module_hamilt_lcao/hamilt_lcaodft/spar_exx.h"

namespace sparse_format
{

    template<typename Tdata> void cal_HR_exx(
            LCAO_Matrix &lm,
            const int &current_spin,
            const double &sparse_thr,
            const int (&nmp)[3],
            const std::vector< std::map <int, std::map < std::pair<int, std::array<int,3>>, RI::Tensor<Tdata> > >>& Hexxs);

}
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_hamilt.hpp"
#endif
#endif
