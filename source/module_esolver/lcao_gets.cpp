#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_esolver/esolver_ks_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
//
#include "module_base/timer.h"
#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_io/write_HS_R.h"
#include "module_parameter/parameter.h"
#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#endif
#include "module_base/formatter.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_hamilt_general/module_ewald/H_Ewald_pw.h"
#include "module_hamilt_general/module_vdw/vdw.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_domain.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/op_exx_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/operator_lcao.h"
#include "module_hamilt_lcao/module_deltaspin/spin_constrain.h"
#include "module_io/read_wfc_nao.h"
#include "module_io/write_elecstat_pot.h"
#include "module_io/write_wfc_nao.h"
#ifdef __EXX
#include "module_io/restart_exx_csr.h"
#endif

namespace ModuleESolver
{

template <>
void ESolver_KS_LCAO<double, double>::get_S(void)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "get_S");
    ModuleBase::WARNING_QUIT("ESolver_KS_LCAO<double,double>::get_S", "not implemented for");
}

template <>
void ESolver_KS_LCAO<std::complex<double>, double>::get_S(void)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "get_S");
    // (1) Find adjacent atoms for each atom.
    double search_radius = -1.0;
    search_radius = atom_arrange::set_sr_NL(GlobalV::ofs_running,
                                                     PARAM.inp.out_level,
                                                     orb_.get_rcutmax_Phi(),
                                                     GlobalC::ucell.infoNL.get_rcutmax_Beta(),
                                                     PARAM.globalv.gamma_only_local);

    atom_arrange::search(PARAM.inp.search_pbc,
                         GlobalV::ofs_running,
                         GlobalC::GridD,
                         GlobalC::ucell,
                         search_radius,
                         PARAM.inp.test_atom_input);

    this->RA.for_2d(this->pv, PARAM.globalv.gamma_only_local, orb_.cutoffs());

    if (this->p_hamilt == nullptr) {
        this->p_hamilt = new hamilt::HamiltLCAO<std::complex<double>, double>(
            &this->pv,
            this->kv,
            *(two_center_bundle_.overlap_orb),
            orb_.cutoffs());
        dynamic_cast<hamilt::OperatorLCAO<std::complex<double>, double>*>(
            this->p_hamilt->ops)
            ->contributeHR();
    }

    // mohan add 2024-06-09
    const std::string fn = PARAM.globalv.global_out_dir + "SR.csr";

    std::cout << " The file is saved in " << fn << std::endl;

    ModuleIO::output_SR(pv, GlobalC::GridD, this->p_hamilt, fn);

    return;
}

template <>
void ESolver_KS_LCAO<std::complex<double>, std::complex<double>>::get_S(void)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "get_S");
    // (1) Find adjacent atoms for each atom.
    double search_radius = -1.0;
    search_radius = atom_arrange::set_sr_NL(GlobalV::ofs_running,
                                                     PARAM.inp.out_level,
                                                     orb_.get_rcutmax_Phi(),
                                                     GlobalC::ucell.infoNL.get_rcutmax_Beta(),
                                                     PARAM.globalv.gamma_only_local);

    atom_arrange::search(PARAM.inp.search_pbc,
                         GlobalV::ofs_running,
                         GlobalC::GridD,
                         GlobalC::ucell,
                         search_radius,
                         PARAM.inp.test_atom_input);

    this->RA.for_2d(this->pv, PARAM.globalv.gamma_only_local, orb_.cutoffs());
    if (this->p_hamilt == nullptr) {
        this->p_hamilt = new hamilt::HamiltLCAO<std::complex<double>,
                                                std::complex<double>>(
            &this->pv,
            this->kv,
            *(two_center_bundle_.overlap_orb),
            orb_.cutoffs()
            );
        dynamic_cast<
            hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>*>(
            this->p_hamilt->ops)
            ->contributeHR();
    }

    // mohan add 2024-06-09
    const std::string fn = PARAM.globalv.global_out_dir + "SR.csr";

    std::cout << " The file is saved in " << fn << std::endl;

    ModuleIO::output_SR(pv, GlobalC::GridD, this->p_hamilt, fn);

    return;
}

template class ESolver_KS_LCAO<double, double>;
template class ESolver_KS_LCAO<std::complex<double>, double>;
template class ESolver_KS_LCAO<std::complex<double>, std::complex<double>>;
} // namespace ModuleESolver
