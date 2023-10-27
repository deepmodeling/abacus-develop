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
    const int nlocal = nw/2;
    ModuleBase::matrix MecMulP, orbMulP;
    MecMulP.create(this->nspin_, nlocal, true);
    orbMulP.create(this->nspin_, nlocal, true);

    for(size_t ik = 0; ik != this->kv_.nks; ++ik)
    {
        dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>*>(this->p_hamilt)
            ->updateSk(ik, LM, 1);

        ModuleBase::ComplexMatrix mud;
        mud.create(this->ParaV->ncol, this->ParaV->nrow, true);

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
                LM->Sloc2.data(),
                &one_int,
                &one_int,
                this->ParaV->desc,
                &zero_float,
                mud.c,
                &one_int,
                &one_int,
                this->ParaV->desc);

        for(size_t i=0; i < nw; ++i)
        {
            const int index = i%2;
            if(!index)
            {
                const int j = i/2;
                const int k1 = 2*j;
                const int k2 = 2*j+1;
                if(this->ParaV->in_this_processor(k1, k1))
                {
                    const int ir = this->ParaV->global2local_row(k1);
                    const int ic = this->ParaV->global2local_col(k1);
                    MecMulP(0, j) += mud(ic, ir).real();
                    MecMulP(3, j) += mud(ic, ir).real();
                }
                if(this->ParaV->in_this_processor(k1, k2))
                {
                    const int ir = this->ParaV->global2local_row(k1);
                    const int ic = this->ParaV->global2local_col(k2);
                    MecMulP(1, j) += mud(ic, ir).real();
                    MecMulP(2, j) += mud(ic, ir).imag();
                }
                if(this->ParaV->in_this_processor(k2, k1))
                {
                    const int ir = this->ParaV->global2local_row(k2);
                    const int ic = this->ParaV->global2local_col(k1);
                    MecMulP(1, j) += mud(ic, ir).real();
                    MecMulP(2, j) -= mud(ic, ir).imag();
                }
                if(this->ParaV->in_this_processor(k2, k2))
                {
                    const int ir = this->ParaV->global2local_row(k2);
                    const int ic = this->ParaV->global2local_col(k2);
                    MecMulP(0, j) += mud(ic, ir).real();
                    MecMulP(3, j) -= mud(ic, ir).real();
                }
            }
        }
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