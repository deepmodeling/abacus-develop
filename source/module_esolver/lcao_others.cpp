#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_esolver/esolver_ks_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
//
#include "module_base/timer.h"
#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_io/berryphase.h"
#include "module_io/get_pchg_lcao.h"
#include "module_io/get_wf_lcao.h"
#include "module_io/to_wannier90_lcao.h"
#include "module_io/to_wannier90_lcao_in_pw.h"
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

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::others(const int istep)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "others");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "others");

    const std::string cal_type = PARAM.inp.calculation;

    if (cal_type == "test_memory")
    {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "testing memory");
        Cal_Test::test_memory(this->pw_rho,
                              this->pw_wfc,
                              this->p_chgmix->get_mixing_mode(),
                              this->p_chgmix->get_mixing_ndim());
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "testing memory");
        return;
    }
    else if (cal_type == "test_neighbour")
    {
        // test_search_neighbor();
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "testing neighbour");
        double search_radius = PARAM.inp.search_radius;
        atom_arrange::search(PARAM.inp.search_pbc,
                             GlobalV::ofs_running,
                             GlobalC::GridD,
                             GlobalC::ucell,
                             search_radius,
                             PARAM.inp.test_atom_input,
                             true);
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "testing neighbour");
        return;
    }

    this->beforesolver(istep);
    // pelec should be initialized before these calculations
    this->pelec->init_scf(istep, this->sf.strucFac, GlobalC::ucell.symm);
    // self consistent calculations for electronic ground state
    if (cal_type == "get_pchg")
    {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "getting partial charge");
        IState_Charge ISC(this->psi, &(this->pv));
        if (PARAM.globalv.gamma_only_local)
        {
            ISC.begin(this->GG,
                      this->pelec->charge->rho,
                      this->pelec->wg,
                      this->pelec->eferm.get_all_ef(),
                      this->pw_rhod->nrxx,
                      this->pw_rhod->nplane,
                      this->pw_rhod->startz_current,
                      this->pw_rhod->nx,
                      this->pw_rhod->ny,
                      this->pw_rhod->nz,
                      this->pw_big->bz,
                      this->pw_big->nbz,
                      PARAM.globalv.gamma_only_local,
                      PARAM.inp.nbands_istate,
                      PARAM.inp.out_pchg,
                      PARAM.inp.nbands,
                      PARAM.inp.nelec,
                      PARAM.inp.nspin,
                      PARAM.globalv.nlocal,
                      PARAM.globalv.global_out_dir,
                      GlobalV::ofs_warning,
                      &GlobalC::ucell,
                      &GlobalC::GridD,
                      this->kv);
        }
        else
        {
            ISC.begin(this->GK,
                      this->pelec->charge->rho,
                      this->pelec->charge->rhog,
                      this->pelec->wg,
                      this->pelec->eferm.get_all_ef(),
                      this->pw_rhod,
                      this->pw_rhod->nrxx,
                      this->pw_rhod->nplane,
                      this->pw_rhod->startz_current,
                      this->pw_rhod->nx,
                      this->pw_rhod->ny,
                      this->pw_rhod->nz,
                      this->pw_big->bz,
                      this->pw_big->nbz,
                      PARAM.globalv.gamma_only_local,
                      PARAM.inp.nbands_istate,
                      PARAM.inp.out_pchg,
                      PARAM.inp.nbands,
                      PARAM.inp.nelec,
                      PARAM.inp.nspin,
                      PARAM.globalv.nlocal,
                      PARAM.globalv.global_out_dir,
                      GlobalV::ofs_warning,
                      &GlobalC::ucell,
                      &GlobalC::GridD,
                      this->kv,
                      PARAM.inp.if_separate_k,
                      &GlobalC::Pgrid,
                      this->pelec->charge->ngmc);
        }
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "getting partial charge");
    }
    else if (cal_type == "get_wf")
    {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "getting wave function");
        IState_Envelope IEP(this->pelec);
        if (PARAM.globalv.gamma_only_local)
        {
            IEP.begin(this->psi,
                      this->pw_rhod,
                      this->pw_wfc,
                      this->pw_big,
                      this->pv,
                      this->GG,
                      PARAM.inp.out_wfc_pw,
                      this->wf.out_wfc_r,
                      this->kv,
                      PARAM.inp.nelec,
                      PARAM.inp.nbands_istate,
                      PARAM.inp.out_wfc_norm,
                      PARAM.inp.out_wfc_re_im,
                      PARAM.inp.nbands,
                      PARAM.inp.nspin,
                      PARAM.globalv.nlocal,
                      PARAM.globalv.global_out_dir);
        }
        else
        {
            IEP.begin(this->psi,
                      this->pw_rhod,
                      this->pw_wfc,
                      this->pw_big,
                      this->pv,
                      this->GK,
                      PARAM.inp.out_wfc_pw,
                      this->wf.out_wfc_r,
                      this->kv,
                      PARAM.inp.nelec,
                      PARAM.inp.nbands_istate,
                      PARAM.inp.out_wfc_norm,
                      PARAM.inp.out_wfc_re_im,
                      PARAM.inp.nbands,
                      PARAM.inp.nspin,
                      PARAM.globalv.nlocal,
                      PARAM.globalv.global_out_dir);
        }
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "getting wave function");
    }
    else
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_LCAO<TK, TR>::others", "CALCULATION type not supported");
    }

    ModuleBase::timer::tick("ESolver_KS_LCAO", "others");
    return;
}

template class ESolver_KS_LCAO<double, double>;
template class ESolver_KS_LCAO<std::complex<double>, double>;
template class ESolver_KS_LCAO<std::complex<double>, std::complex<double>>;
} // namespace ModuleESolver
