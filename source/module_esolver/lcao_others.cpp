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

    // 1. prepare HS matrices, prepare grid integral
    // (1) Find adjacent atoms for each atom.
    double search_radius = atom_arrange::set_sr_NL(GlobalV::ofs_running,
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

    // (3) Periodic condition search for each grid.
    double dr_uniform = 0.001;
    std::vector<double> rcuts;
    std::vector<std::vector<double>> psi_u;
    std::vector<std::vector<double>> dpsi_u;
    std::vector<std::vector<double>> d2psi_u;

    Gint_Tools::init_orb(dr_uniform, rcuts, GlobalC::ucell, orb_, psi_u, dpsi_u, d2psi_u);

    this->GridT.set_pbc_grid(this->pw_rho->nx,
                             this->pw_rho->ny,
                             this->pw_rho->nz,
                             this->pw_big->bx,
                             this->pw_big->by,
                             this->pw_big->bz,
                             this->pw_big->nbx,
                             this->pw_big->nby,
                             this->pw_big->nbz,
                             this->pw_big->nbxx,
                             this->pw_big->nbzp_start,
                             this->pw_big->nbzp,
                             this->pw_rho->ny,
                             this->pw_rho->nplane,
                             this->pw_rho->startz_current,
                             GlobalC::ucell,
                             GlobalC::GridD,
                             dr_uniform,
                             rcuts,
                             psi_u,
                             dpsi_u,
                             d2psi_u,
                             PARAM.inp.nstream);
    psi_u.clear();
    psi_u.shrink_to_fit();
    dpsi_u.clear();
    dpsi_u.shrink_to_fit();
    d2psi_u.clear();
    d2psi_u.shrink_to_fit();

    // (2)For each atom, calculate the adjacent atoms in different cells
    // and allocate the space for H(R) and S(R).
    // If k point is used here, allocate HlocR after atom_arrange.
    this->RA.for_2d(this->pv, PARAM.globalv.gamma_only_local, orb_.cutoffs());

    // 2. density matrix extrapolation

    // set the augmented orbitals index.
    // after ParaO and GridT,
    // this information is used to calculate
    // the force.

    // init psi
    if (this->psi == nullptr)
    {
        int nsk = 0;
        int ncol = 0;
        if (PARAM.globalv.gamma_only_local)
        {
            nsk = PARAM.inp.nspin;
            ncol = this->pv.ncol_bands;
            if (PARAM.inp.ks_solver == "genelpa" || PARAM.inp.ks_solver == "elpa" || PARAM.inp.ks_solver == "lapack"
                || PARAM.inp.ks_solver == "pexsi" || PARAM.inp.ks_solver == "cusolver"
                || PARAM.inp.ks_solver == "cusolvermp")
            {
                ncol = this->pv.ncol;
            }
        }
        else
        {
            nsk = this->kv.get_nks();
#ifdef __MPI
            ncol = this->pv.ncol_bands;
#else
            ncol = PARAM.inp.nbands;
#endif
        }
        this->psi = new psi::Psi<TK>(nsk, ncol, this->pv.nrow, nullptr);
    }

    // init wfc from file
    if (istep == 0 && PARAM.inp.init_wfc == "file")
    {
        if (!ModuleIO::read_wfc_nao(PARAM.globalv.global_readin_dir, this->pv, *(this->psi), this->pelec))
        {
            ModuleBase::WARNING_QUIT("ESolver_KS_LCAO<TK, TR>::beforesolver", "read wfc nao failed");
        }
    }

    // prepare grid in Gint
    LCAO_domain::grid_prepare(this->GridT, this->GG, this->GK, orb_, *this->pw_rho, *this->pw_big);

    // init Hamiltonian
    if (this->p_hamilt != nullptr)
    {
        delete this->p_hamilt;
        this->p_hamilt = nullptr;
    }
    if (this->p_hamilt == nullptr)
    {
        elecstate::DensityMatrix<TK, double>* DM = dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM();
        this->p_hamilt = new hamilt::HamiltLCAO<TK, TR>(
            PARAM.globalv.gamma_only_local ? &(this->GG) : nullptr,
            PARAM.globalv.gamma_only_local ? nullptr : &(this->GK),
            &this->pv,
            this->pelec->pot,
            this->kv,
            two_center_bundle_,
            orb_,
            DM
#ifdef __EXX
            ,
            istep,
            GlobalC::exx_info.info_ri.real_number ? &this->exd->two_level_step : &this->exc->two_level_step,
            GlobalC::exx_info.info_ri.real_number ? &exx_lri_double->Hexxs : nullptr,
            GlobalC::exx_info.info_ri.real_number ? nullptr : &exx_lri_complex->Hexxs
#endif
        );
    }

#ifdef __DEEPKS
    // for each ionic step, the overlap <psi|alpha> must be rebuilt
    // since it depends on ionic positions
    if (PARAM.globalv.deepks_setorb)
    {
        const Parallel_Orbitals* pv = &this->pv;
        // build and save <psi(0)|alpha(R)> at beginning
        GlobalC::ld.build_psialpha(PARAM.inp.cal_force,
                                   GlobalC::ucell,
                                   orb_,
                                   GlobalC::GridD,
                                   *(two_center_bundle_.overlap_orb_alpha));

        if (PARAM.inp.deepks_out_unittest)
        {
            GlobalC::ld.check_psialpha(PARAM.inp.cal_force, GlobalC::ucell, orb_, GlobalC::GridD);
        }
    }
#endif
    if (PARAM.inp.sc_mag_switch)
    {
        spinconstrain::SpinConstrain<TK>& sc = spinconstrain::SpinConstrain<TK>::getScInstance();
        sc.init_sc(PARAM.inp.sc_thr,
                   PARAM.inp.nsc,
                   PARAM.inp.nsc_min,
                   PARAM.inp.alpha_trial,
                   PARAM.inp.sccut,
                   PARAM.inp.sc_drop_thr,
                   GlobalC::ucell,
                   &(this->pv),
                   PARAM.inp.nspin,
                   this->kv,
                   PARAM.inp.ks_solver,
                   this->p_hamilt,
                   this->psi,
                   this->pelec);
    }
    //=========================================================
    // cal_ux should be called before init_scf because
    // the direction of ux is used in noncoline_rho
    //=========================================================
    if (PARAM.inp.nspin == 4)
    {
        GlobalC::ucell.cal_ux();
    }

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
