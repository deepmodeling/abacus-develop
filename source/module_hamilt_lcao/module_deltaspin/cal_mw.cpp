#include <iostream>

#include "module_base/matrix.h"
#include "module_base/name_angular.h"
#include "module_base/scalapack_connector.h"
#include "module_base/tool_title.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "spin_constrain.h"

template <>
ModuleBase::matrix SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::cal_MW_k(
    LCAO_Matrix* LM,
    const std::vector<std::vector<std::complex<double>>>& dm)
{
    ModuleBase::TITLE("module_deltaspin", "cal_MW_k");
    int nw = this->get_nw();
    const int nlocal = (this->nspin_ == 4) ? nw / 2 : nw;
    ModuleBase::matrix MecMulP(this->nspin_, nlocal, true), orbMulP(this->nspin_, nlocal, true);
    for(size_t ik = 0; ik != this->kv_.nks; ++ik)
    {
        if (this->nspin_ == 4)
        {
            dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>*>(this->p_hamilt)->updateSk(ik, LM, 1);
        }
        else
        {
            dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, double>*>(this->p_hamilt)->updateSk(ik, LM, 1);
        }
        // calculat weight function
        std::vector<std::complex<double>> weight(LM->ParaV->nloc);
        std::fill(weight.begin(), weight.end(), std::complex<double>(0, 0));
        this->cal_weight_func(&weight[0],
                              LM->Sloc2,
                              ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(),
                              this->kv_.isk[ik]);
        ModuleBase::ComplexMatrix mud(this->ParaV->ncol, this->ParaV->nrow, true);
#ifdef __MPI
        const char T_char = 'T';
        const char N_char = 'N';
        const int one_int = 1;
        const std::complex<double> one_float = {1.0, 0.0}, zero_float = {0.0, 0.0};
        pzgemm_(&T_char,
                &T_char,
                &nw,
                &nw,
                &nw,
                &one_float,
                dm[ik].data(),
                &one_int,
                &one_int,
                this->ParaV->desc,
                weight.data(),
                &one_int,
                &one_int,
                this->ParaV->desc,
                &zero_float,
                mud.c,
                &one_int,
                &one_int,
                this->ParaV->desc);
        this->collect_MW(MecMulP, mud, nw, this->kv_.isk[ik]);
#endif
    }
#ifdef __MPI
    MPI_Allreduce(MecMulP.c, orbMulP.c, this->nspin_*nlocal, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif 

    return orbMulP;
}

template <>
void SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::cal_MW(const int& step,
                                                                  LCAO_Matrix* LM,
                                                                  bool print)
{
    ModuleBase::TITLE("module_deltaspin", "cal_MW");
    const std::vector<std::vector<std::complex<double>>>& dm
        = dynamic_cast<const elecstate::ElecStateLCAO<std::complex<double>>*>(this->pelec)->get_DM()->get_DMK_vector();
    this->calculate_MW(this->convert(this->cal_MW_k(LM, dm)));
    this->print_Mi(print);
}