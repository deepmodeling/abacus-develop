#include "source_estate/elecstate_lcao.h"
#include "source_estate/cal_dm.h"
#include "source_base/timer.h"
#include "source_estate/module_dm/cal_dm_psi.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_lcao/module_deltaspin/spin_constrain.h"
#include "source_io/module_parameter/parameter.h"

#include "source_lcao/module_gint/gint_interface.h"

#include <vector>

namespace elecstate
{


template <>
double ElecStateLCAO<double>::get_spin_constrain_energy()
{
    spinconstrain::SpinConstrain<double>& sc = spinconstrain::SpinConstrain<double>::getScInstance();
    return sc.cal_escon();
}

template <>
double ElecStateLCAO<std::complex<double>>::get_spin_constrain_energy()
{
    spinconstrain::SpinConstrain<std::complex<double>>& sc
        = spinconstrain::SpinConstrain<std::complex<double>>::getScInstance();
    return sc.cal_escon();
}

template <>
void ElecStateLCAO<double>::dm2rho(std::vector<double*> pexsi_DM, 
		std::vector<double*> pexsi_EDM,
		DensityMatrix<double, double>* dm)
{
    ModuleBase::timer::start("ElecStateLCAO", "dm2rho");

    int nspin = PARAM.inp.nspin;
    if (PARAM.inp.nspin == 4)
    {
        nspin = 1;
    }

#ifdef __PEXSI
    dm->pexsi_EDM = pexsi_EDM;
#endif

    for (int is = 0; is < nspin; is++)
    {
        dm->set_DMK_pointer(is, pexsi_DM[is]);
    }
    dm->cal_DMR();

    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        ModuleBase::GlobalFunc::ZEROS(this->charge->rho[is],
                                      this->charge->nrxx); // mohan 2009-11-10
    }

    ModuleBase::GlobalFunc::NOTE("Calculate the charge on real space grid!");
    ModuleGint::cal_gint_rho(dm->get_DMR_vector(), PARAM.inp.nspin, this->charge->rho);
    if (XC_Functional::get_ked_flag())
    {
        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(this->charge->kin_r[0], this->charge->nrxx);
        }
        ModuleGint::cal_gint_tau(dm->get_DMR_vector(), PARAM.inp.nspin, this->charge->kin_r);
    }

    this->charge->renormalize_rho();

    ModuleBase::timer::end("ElecStateLCAO", "dm2rho");
    return;
}

template <>
void ElecStateLCAO<std::complex<double>>::dm2rho(std::vector<std::complex<double>*> pexsi_DM,
		std::vector<std::complex<double>*> pexsi_EDM,
		DensityMatrix<std::complex<double>, double>* dm)
{
    ModuleBase::timer::start("ElecStateLCAO", "dm2rho");

    const int dmk_size = dm->get_DMK_size();
    if (static_cast<int>(pexsi_DM.size()) < dmk_size)
    {
        ModuleBase::WARNING_QUIT("ElecStateLCAO", "PEXSI multi-k density matrix buffer is smaller than DMK storage");
    }

    const int local_matrix_size = dm->get_DMK_nrow() * dm->get_DMK_ncol();
    const double pexsi_spin_weight = PARAM.inp.nspin == 1 ? 0.5 : 1.0;
    for (int ik = 0; ik < dmk_size; ++ik)
    {
        const double k_weight = ((this->klist != nullptr && ik < static_cast<int>(this->klist->wk.size()))
                                     ? this->klist->wk[ik]
                                     : 1.0)
                                * pexsi_spin_weight;
        for (int i = 0; i < local_matrix_size; ++i)
        {
            pexsi_DM[ik][i] *= k_weight;
            pexsi_EDM[ik][i] *= k_weight;
        }
        dm->set_DMK_pointer(ik, pexsi_DM[ik]);
    }

#ifdef __PEXSI
    dm->pexsi_EDM = pexsi_EDM;
#endif

    dm->cal_DMR();

    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        ModuleBase::GlobalFunc::ZEROS(this->charge->rho[is], this->charge->nrxx);
    }

    ModuleBase::GlobalFunc::NOTE("Calculate the charge on real space grid!");
    ModuleGint::cal_gint_rho(dm->get_DMR_vector(), PARAM.inp.nspin, this->charge->rho);
    if (XC_Functional::get_ked_flag())
    {
        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
            ModuleBase::GlobalFunc::ZEROS(this->charge->kin_r[0], this->charge->nrxx);
        }
        ModuleGint::cal_gint_tau(dm->get_DMR_vector(), PARAM.inp.nspin, this->charge->kin_r);
    }

    this->charge->renormalize_rho();

    ModuleBase::timer::end("ElecStateLCAO", "dm2rho");
    return;
}


template class ElecStateLCAO<double>;               // Gamma_only case
template class ElecStateLCAO<std::complex<double>>; // multi-k case

} // namespace elecstate
