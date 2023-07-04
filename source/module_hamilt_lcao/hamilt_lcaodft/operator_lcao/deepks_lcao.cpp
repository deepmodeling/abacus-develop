#include "deepks_lcao.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#endif
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"

namespace hamilt
{

template class DeePKS<OperatorLCAO<double>>;

template class DeePKS<OperatorLCAO<std::complex<double>>>;

template<>
void DeePKS<OperatorLCAO<double>>::contributeHR()
{
    ModuleBase::TITLE("DeePKS", "contributeHR");
#ifdef __DEEPKS
    if(GlobalC::ld.get_hr_cal())
    {
        ModuleBase::timer::tick("DeePKS", "contributeHR");
        const Parallel_Orbitals* pv = this->LM->ParaV;
        GlobalC::ld.cal_projected_DM(this->loc->dm_gamma[0],
            GlobalC::ucell,
            GlobalC::ORB,
            GlobalC::GridD,
            pv->global2local_row,
            pv->global2local_col);
        GlobalC::ld.cal_descriptor();        
        GlobalC::ld.cal_gedm(GlobalC::ucell.nat);
        GlobalC::ld.add_v_delta(GlobalC::ucell,
            GlobalC::ORB,
            GlobalC::GridD,
            pv->global2local_row,
            pv->global2local_col,
            pv->nrow,
            pv->ncol);

        GlobalC::ld.set_hr_cal(false);

        ModuleBase::timer::tick("DeePKS", "contributeHR");
    }
#endif
}

template<>
void DeePKS<OperatorLCAO<std::complex<double>>>::contributeHR()
{
#ifdef __DEEPKS
    ModuleBase::TITLE("DeePKS", "contributeHR");
    // if DM_K changed, HR of DeePKS need to refresh.
    // the judgement is based on the status of HR in GlobalC::ld
    // this operator should be informed that DM_K has changed and HR need to recalculate.
    if(GlobalC::ld.get_hr_cal())
    {
        ModuleBase::timer::tick("DeePKS", "contributeHR");

        GlobalC::ld.cal_projected_DM_k(this->loc->dm_k,
                                        GlobalC::ucell,
                                        GlobalC::ORB,
                                        GlobalC::GridD,
                                        this->LM->ParaV->global2local_row,
                                        this->LM->ParaV->global2local_col,
                                        this->nks,
                                        this->kvec_d);
        GlobalC::ld.cal_descriptor();
        // calculate dE/dD
        GlobalC::ld.cal_gedm(GlobalC::ucell.nat);

        // calculate H_V_deltaR from saved <alpha(0)|psi(R)>
        GlobalC::ld
            .add_v_delta_k(GlobalC::ucell, GlobalC::ORB, GlobalC::GridD, this->LM->ParaV->global2local_row, this->LM->ParaV->global2local_col, this->LM->ParaV->nnr);
        
        GlobalC::ld.set_hr_cal(false);
        
        ModuleBase::timer::tick("DeePKS", "contributeHR");
    } 

#endif
}

template<>
void DeePKS<OperatorLCAO<double>>::contributeHk(int ik)
{
#ifdef __DEEPKS	//caoyu add 2021-07-26 for DeePKS

    ModuleBase::TITLE("DeePKS", "contributeHk");
    ModuleBase::timer::tick("DeePKS", "contributeHk");
    
    
    for(int iic=0;iic<this->LM->ParaV->nloc;iic++)
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