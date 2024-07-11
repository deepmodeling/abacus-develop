#ifndef LCAO_MATRIX_H
#define LCAO_MATRIX_H

#include "module_base/complexmatrix.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/vector3.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_HS_arrays.hpp"

// add by jingan for map<> in 2021-12-2, will be deleted in the future
#include "module_base/abfs-vector3_order.h"
#ifdef __EXX
#include <RI/global/Tensor.h>
#endif

class LCAO_Matrix {
  public:
   LCAO_Matrix(){};
  ~LCAO_Matrix(){};

    Parallel_Orbitals* ParaV;

#ifdef __EXX
    using TAC = std::pair<int, std::array<int, 3>>;
    std::vector<std::map<int, std::map<TAC, RI::Tensor<double>>>>* Hexxd;
    std::vector<std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>*
        Hexxc;
    /// @brief Hexxk for all k-points, only for the 1st scf loop ofrestart load
    std::vector<std::vector<double>> Hexxd_k_load;
    std::vector<std::vector<std::complex<double>>> Hexxc_k_load;
#endif
};

#endif
