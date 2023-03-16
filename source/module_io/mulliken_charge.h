#ifndef MULLIKEN_CHARGE_H
#define MULLIKEN_CHARGE_H

#include "module_base/matrix.h"
#include "module_base/complexmatrix.h"
#ifdef __LCAO
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_hamilt.h"
#endif

// by qifeng, refactor by jiyy 2023-02-25
class Mulliken_Charge
{
	public:

	void out_mulliken(LCAO_Hamilt &uhm, Local_Orbital_Charge &loc);

	private:
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
        LCAO_Hamilt &uhm
    );

    std::vector<std::vector<std::vector<double>>> convert(const ModuleBase::matrix &orbMulP);

    double output_cut(const double& result);
};
#endif