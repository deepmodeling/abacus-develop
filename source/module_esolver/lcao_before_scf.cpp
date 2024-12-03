#include "module_elecstate/cal_ux.h"
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
#include "module_io/io_npz.h"
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
#include "module_io/cube_io.h"
#include "module_io/read_wfc_nao.h"
#include "module_io/write_elecstat_pot.h"
#include "module_io/write_wfc_nao.h"
#ifdef __EXX
#include "module_io/restart_exx_csr.h"
#endif

namespace ModuleESolver
{

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::before_scf(UnitCell& ucell, const int istep)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "before_scf");

    //! 1) call before_scf() of ESolver_FP
    ESolver_FP::before_scf(ucell, istep);

    if (ucell.ionic_position_updated)
    {
        this->CE.update_all_dis(ucell);
        this->CE.extrapolate_charge(
#ifdef __MPI
            &(GlobalC::Pgrid),
#endif
            ucell,
            this->pelec->charge,
            &(this->sf),
            GlobalV::ofs_running,
            GlobalV::ofs_warning);
    }

    //----------------------------------------------------------
    // about vdw, jiyy add vdwd3 and linpz add vdwd2
    //----------------------------------------------------------
    auto vdw_solver = vdw::make_vdw(ucell, PARAM.inp, &(GlobalV::ofs_running));
    if (vdw_solver != nullptr)
    {
        this->pelec->f_en.evdw = vdw_solver->get_energy();
    }

    // 1. prepare HS matrices, prepare grid integral
    // (1) Find adjacent atoms for each atom.
    double search_radius = atom_arrange::set_sr_NL(GlobalV::ofs_running,
                                                   PARAM.inp.out_level,
                                                   orb_.get_rcutmax_Phi(),
                                                   ucell.infoNL.get_rcutmax_Beta(),
                                                   PARAM.globalv.gamma_only_local);

    atom_arrange::search(PARAM.inp.search_pbc,
                         GlobalV::ofs_running,
                         GlobalC::GridD,
                         ucell,
                         search_radius,
                         PARAM.inp.test_atom_input);

    // (3) Periodic condition search for each grid.
    double dr_uniform = 0.001;
    std::vector<double> rcuts;
    std::vector<std::vector<double>> psi_u;
    std::vector<std::vector<double>> dpsi_u;
    std::vector<std::vector<double>> d2psi_u;

    Gint_Tools::init_orb(dr_uniform, rcuts, ucell, orb_, psi_u, dpsi_u, d2psi_u);

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
                             ucell,
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
    this->RA.for_2d(ucell, GlobalC::GridD, this->pv, PARAM.globalv.gamma_only_local, orb_.cutoffs());

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
    LCAO_domain::grid_prepare(this->GridT, this->GG, this->GK, ucell, orb_, *this->pw_rho, *this->pw_big);

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
            ucell,
            GlobalC::GridD,
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
                                   ucell,
                                   orb_,
                                   GlobalC::GridD,
                                   *(two_center_bundle_.overlap_orb_alpha));

        if (PARAM.inp.deepks_out_unittest)
        {
            GlobalC::ld.check_psialpha(PARAM.inp.cal_force, ucell, orb_, GlobalC::GridD);
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
                   ucell,
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
    elecstate::cal_ux(ucell);

    // Peize Lin add 2016-12-03
#ifdef __EXX // set xc type before the first cal of xc in pelec->init_scf
    if (PARAM.inp.calculation != "nscf")
    {
        if (GlobalC::exx_info.info_ri.real_number)
        {
            this->exd->exx_beforescf(istep, this->kv, *this->p_chgmix, ucell, orb_);
        }
        else
        {
            this->exc->exx_beforescf(istep, this->kv, *this->p_chgmix, ucell, orb_);
        }
    }
#endif // __EXX

    this->pelec->init_scf(istep, this->sf.strucFac, this->ppcell.numeric, ucell.symm);

    //! output the initial charge density
    if (PARAM.inp.out_chg[0] == 2)
    {
        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
            std::stringstream ss;
            ss << PARAM.globalv.global_out_dir << "SPIN" << is + 1 << "_CHG_INI.cube";
            ModuleIO::write_vdata_palgrid(GlobalC::Pgrid,
                                          this->pelec->charge->rho[is],
                                          is,
                                          PARAM.inp.nspin,
                                          istep,
                                          ss.str(),
                                          this->pelec->eferm.ef,
                                          &(ucell));
        }
    }

    //! output total local potential of the initial charge density
    if (PARAM.inp.out_pot == 3)
    {
        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
            std::stringstream ss;
            ss << PARAM.globalv.global_out_dir << "SPIN" << is + 1 << "_POT_INI.cube";
            ModuleIO::write_vdata_palgrid(GlobalC::Pgrid,
                                          this->pelec->pot->get_effective_v(is),
                                          is,
                                          PARAM.inp.nspin,
                                          istep,
                                          ss.str(),
                                          0.0, // efermi
                                          &(ucell),
                                          11, // precsion
                                          0); // out_fermi
        }
    }

    // initalize DMR
    // DMR should be same size with Hamiltonian(R)
    dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)
        ->get_DM()
        ->init_DMR(*(dynamic_cast<hamilt::HamiltLCAO<TK, TR>*>(this->p_hamilt)->getHR()));
    // two cases are considered:
    // 1. DMK in DensityMatrix is not empty (istep > 0), then DMR is initialized by DMK
    // 2. DMK in DensityMatrix is empty (istep == 0), then DMR is initialized by zeros
    if (istep > 0)
    {
        dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM()->cal_DMR();
    }

    if (PARAM.inp.dm_to_rho)
    {
        std::string zipname = "output_DM0.npz";
        elecstate::DensityMatrix<TK, double>* dm
            = dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM();
        ModuleIO::read_mat_npz(&(this->pv), ucell, zipname, *(dm->get_DMR_pointer(1)));
        if (PARAM.inp.nspin == 2)
        {
            zipname = "output_DM1.npz";
            ModuleIO::read_mat_npz(&(this->pv), ucell, zipname, *(dm->get_DMR_pointer(2)));
        }

        this->pelec->calculate_weights();
        this->pelec->psiToRho(*this->psi);

        int nspin0 = PARAM.inp.nspin == 2 ? 2 : 1;
        for (int is = 0; is < nspin0; is++)
        {
            std::string fn = PARAM.globalv.global_out_dir + "/SPIN" + std::to_string(is + 1) + "_CHG.cube";
            ModuleIO::write_vdata_palgrid(GlobalC::Pgrid,
                                          this->pelec->charge->rho[is],
                                          is,
                                          PARAM.inp.nspin,
                                          istep,
                                          fn,
                                          this->pelec->eferm.get_efval(is),
                                          &(ucell),
                                          3,
                                          1);
        }

        return;
    }

    // the electron charge density should be symmetrized,
    // here is the initialization
    Symmetry_rho srho;
    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        srho.begin(is, *(this->pelec->charge), this->pw_rho, ucell.symm);
    }

    // 1. calculate ewald energy.
    // mohan update 2021-02-25
    if (!PARAM.inp.test_skip_ewald)
    {
        this->pelec->f_en.ewald_energy = H_Ewald_pw::compute_ewald(ucell, this->pw_rho, this->sf.strucFac);
    }

    this->p_hamilt->non_first_scf = istep;

    // update in ion-step
    if (PARAM.inp.rdmft == true)
    {
        // necessary operation of these parameters have be done with p_esolver->Init() in source/driver_run.cpp
        rdmft_solver.update_ion(ucell,
                                *(this->pw_rho),
                                this->ppcell.vloc,
                                this->sf.strucFac); // add by jghan, 2024-03-16/2024-10-08
    }

    return;
}

template class ESolver_KS_LCAO<double, double>;
template class ESolver_KS_LCAO<std::complex<double>, double>;
template class ESolver_KS_LCAO<std::complex<double>, std::complex<double>>;
} // namespace ModuleESolver
