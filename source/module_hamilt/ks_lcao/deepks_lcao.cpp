#include "deepks_lcao.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#ifdef __DEEPKS
#include "module_deepks/LCAO_deepks.h"
#endif
#include "src_pw/global.h"

namespace hamilt
{

template class DeePKS<OperatorLCAO<double>>;

template class DeePKS<OperatorLCAO<std::complex<double>>>;

template<>
void DeePKS<OperatorLCAO<double>>::contributeHR()
{

}

template<>
void DeePKS<OperatorLCAO<std::complex<double>>>::contributeHR()
{
#ifdef __DEEPKS
    ModuleBase::TITLE("DeePKS", "contributeHR");
    ModuleBase::timer::tick("DeePKS", "contributeHR");

    GlobalC::ld.cal_projected_DM_k(this->loc->dm_k,
                                    GlobalC::ucell,
                                    GlobalC::ORB,
                                    GlobalC::GridD,
                                    this->LM->ParaV->trace_loc_row,
                                    this->LM->ParaV->trace_loc_col,
                                    GlobalC::kv.nks,
                                    GlobalC::kv.kvec_d);
    GlobalC::ld.cal_descriptor();
    // calculate dE/dD
    GlobalC::ld.cal_gedm(GlobalC::ucell.nat);

    // calculate H_V_deltaR from saved <alpha(0)|psi(R)>
    GlobalC::ld
        .add_v_delta_k(GlobalC::ucell, GlobalC::ORB, GlobalC::GridD, this->LM->ParaV->trace_loc_row, this->LM->ParaV->trace_loc_col, this->LM->ParaV->nnr);
    
    ModuleBase::timer::tick("DeePKS", "contributeHR");

#endif
}

template<>
void DeePKS<OperatorLCAO<double>>::contributeHk(int ik)
{
#ifdef __DEEPKS	//caoyu add 2021-07-26 for DeePKS

    ModuleBase::TITLE("DeePKS", "contributeHk");
    ModuleBase::timer::tick("DeePKS", "contributeHk");
    
    const Parallel_Orbitals* pv = this->LM->ParaV;
    GlobalC::ld.cal_projected_DM(this->loc->dm_gamma[0],
        GlobalC::ucell,
        GlobalC::ORB,
        GlobalC::GridD,
        pv->trace_loc_row,
        pv->trace_loc_col);
    GlobalC::ld.cal_descriptor();        
    GlobalC::ld.cal_gedm(GlobalC::ucell.nat);
    GlobalC::ld.add_v_delta(GlobalC::ucell,
        GlobalC::ORB,
        GlobalC::GridD,
        pv->trace_loc_row,
        pv->trace_loc_col,
        pv->nrow,
        pv->ncol);
    for(int iic=0;iic<pv->nloc;iic++)
    {
        this->LM->Hloc[iic] += GlobalC::ld.H_V_delta[iic];
    }
	
    ModuleBase::timer::tick("DeePKS", "contributeHk");

#endif
}

template<>
void DeePKS<OperatorLCAO<std::complex<double>>>::contributeHk(int ik)
{
    //has been done in folding_fixedH()
}

}