#ifndef MULLIKEN_CHARGE_H
#define MULLIKEN_CHARGE_H

#include "module_base/matrix.h"
#include "module_base/complexmatrix.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_hamilt.h"
#include "module_cell/klist.h"

// by qifeng, refactor by jiyy 2023-02-25
// convert to namespace, liuyu 2023-04-18
namespace ModuleIO
{
	void out_mulliken(const int& step, LCAO_Hamilt &uhm, Local_Orbital_Charge &loc, const K_Vectors& kv);

    /* 
    1. cal_mulliken:    for gamma-only
    2. cal_mulliken_k:  for multi-k

        for `nspin=1` and `nspin=2`:
            return  ModuleBase::matrix with shape (GlobalV::NSPIN, GlobalV::NLOCAL)
        for `nspin=4`:
            return  ModuleBase::matrix with shape (GlobalV::NSPIN, GlobalV::NLOCAL/2)
    */

    ModuleBase::matrix cal_mulliken(const std::vector<ModuleBase::matrix> &dm,
        LCAO_Hamilt &uhm
    );

    ModuleBase::matrix cal_mulliken_k(const std::vector<ModuleBase::ComplexMatrix> &dm,
        LCAO_Hamilt &uhm, const K_Vectors& kv
    );

    std::vector<std::vector<std::vector<double>>> convert(const ModuleBase::matrix &orbMulP);

    inline double output_cut(const double& result)
    {
        if(std::abs(result) < 1e-6)
        {
            return 0.0;
        }
        return result;
    }
}
#endif