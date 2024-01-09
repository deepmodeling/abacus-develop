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
void ESolver_KS_LCAO<TK, TR>::beforesolver(const int istep)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "beforesolver");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "beforesolver");

    // 1. prepare HS matrices, prepare grid integral
    this->set_matrix_grid(this->RA);

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
            nsk = GlobalV::NSPIN;
            ncol = this->pv.ncol_bands;
            if (GlobalV::KS_SOLVER == "genelpa"
                || GlobalV::KS_SOLVER == "elpa"
                || GlobalV::KS_SOLVER == "lapack"
                || GlobalV::KS_SOLVER == "pexsi"
                || GlobalV::KS_SOLVER == "cusolver"
                || GlobalV::KS_SOLVER == "cusolvermp") {
                ncol = this->pv.ncol;
            }
        }
        else
        {
            nsk = this->kv.get_nks();
#ifdef __MPI
            ncol = this->pv.ncol_bands;
#else
            ncol = GlobalV::NBANDS;
#endif
        }
        this->psi = new psi::Psi<TK>(nsk, ncol, this->pv.nrow, nullptr);
    }

    // init wfc from file
    if(istep == 0 && PARAM.inp.init_wfc == "file")
    {
        if (! ModuleIO::read_wfc_nao(GlobalV::global_readin_dir, this->pv, *(this->psi), this->pelec))
        {
            ModuleBase::WARNING_QUIT("ESolver_KS_LCAO<TK, TR>::beforesolver",
                                     "read wfc nao failed");
        }
    }

    // prepare grid in Gint
    LCAO_domain::grid_prepare(this->GridT, this->GG, this->GK, *this->pw_rho, *this->pw_big);

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
            DM
#ifdef __EXX
            , GlobalC::exx_info.info_ri.real_number ? &this->exd->two_level_step : &this->exc->two_level_step
            , GlobalC::exx_info.info_ri.real_number ? &exx_lri_double->Hexxs : nullptr
            , GlobalC::exx_info.info_ri.real_number ? nullptr : &exx_lri_complex->Hexxs
#endif
        );
    }

#ifdef __DEEPKS
    // for each ionic step, the overlap <psi|alpha> must be rebuilt
    // since it depends on ionic positions
    if (GlobalV::deepks_setorb)
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
        SpinConstrain<TK, base_device::DEVICE_CPU>& sc = SpinConstrain<TK, base_device::DEVICE_CPU>::getScInstance();
        sc.init_sc(PARAM.inp.sc_thr,
                   PARAM.inp.nsc,
                   PARAM.inp.nsc_min,
                   PARAM.inp.alpha_trial,
                   PARAM.inp.sccut,
                   PARAM.inp.sc_mag_switch,
                   GlobalC::ucell,
                   PARAM.inp.sc_file,
                   GlobalV::NPOL,
                   &(this->pv),
                   GlobalV::NSPIN,
                   this->kv,
                   GlobalV::KS_SOLVER,
                   this->p_hamilt,
                   this->psi,
                   this->pelec);
    }
    //=========================================================
    // cal_ux should be called before init_scf because
    // the direction of ux is used in noncoline_rho
    //=========================================================
    if (GlobalV::NSPIN == 4 && GlobalV::DOMAG)
    {
        GlobalC::ucell.cal_ux();
    }
    ModuleBase::timer::tick("ESolver_KS_LCAO", "beforesolver");
}

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::before_scf(const int istep)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "before_scf");

    if (GlobalC::ucell.cell_parameter_updated)
    {
        this->init_after_vc(PARAM.inp, GlobalC::ucell);
    }
    if (GlobalC::ucell.ionic_position_updated)
    {
        this->CE.update_all_dis(GlobalC::ucell);
        this->CE.extrapolate_charge(
#ifdef __MPI
            &(GlobalC::Pgrid),
#endif
            GlobalC::ucell,
            this->pelec->charge,
            &(this->sf),
            GlobalV::ofs_running,
            GlobalV::ofs_warning);
    }

    //----------------------------------------------------------
    // about vdw, jiyy add vdwd3 and linpz add vdwd2
    //----------------------------------------------------------
    auto vdw_solver = vdw::make_vdw(GlobalC::ucell, PARAM.inp);
    if (vdw_solver != nullptr)
    {
        this->pelec->f_en.evdw = vdw_solver->get_energy();
    }

    this->beforesolver(istep);

    // Peize Lin add 2016-12-03
#ifdef __EXX // set xc type before the first cal of xc in pelec->init_scf
    if (GlobalC::exx_info.info_ri.real_number)
    {
        this->exd->exx_beforescf(this->kv, *this->p_chgmix, GlobalC::ucell, this->pv);
    }
    else
    {
        this->exc->exx_beforescf(this->kv, *this->p_chgmix, GlobalC::ucell, this->pv);
    }
#endif // __EXX

    this->pelec->init_scf(istep, this->sf.strucFac);

    //! output the initial charge density
    if (PARAM.inp.out_chg[0] == 2)
    {
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            std::stringstream ss;
            ss << GlobalV::global_out_dir << "SPIN" << is + 1 << "_CHG_INI.cube";
            ModuleIO::write_cube(
#ifdef __MPI
                this->pw_big->bz, // bz first, then nbz
                this->pw_big->nbz,
                this->pw_rhod->nplane,
                this->pw_rhod->startz_current,
#endif
                this->pelec->charge->rho[is],
                is,
                GlobalV::NSPIN,
                istep,
                ss.str(),
                this->pw_rhod->nx,
                this->pw_rhod->ny,
                this->pw_rhod->nz,
                this->pelec->eferm.ef,
                &(GlobalC::ucell));
        }
    }

    //! output total local potential of the initial charge density
    if (PARAM.inp.out_pot == 3)
    {
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            std::stringstream ss;
            ss << GlobalV::global_out_dir << "SPIN" << is + 1 << "_POT_INI.cube";
            ModuleIO::write_cube(
#ifdef __MPI
                this->pw_big->bz,
                this->pw_big->nbz,
                this->pw_rhod->nplane,
                this->pw_rhod->startz_current,
#endif
                this->pelec->pot->get_effective_v(is),
                is,
                GlobalV::NSPIN,
                istep,
                ss.str(),
                this->pw_rhod->nx,
                this->pw_rhod->ny,
                this->pw_rhod->nz,
                0.0, // efermi
                &(GlobalC::ucell),
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
    if(istep > 0)
    {
        dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec)
            ->get_DM()
            ->cal_DMR();
    }

    if (PARAM.inp.dm_to_rho)
    {
        std::string zipname = "output_DM0.npz";
        elecstate::DensityMatrix<TK, double>* dm
            = dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM();
        this->read_mat_npz(zipname, *(dm->get_DMR_pointer(1)));
        if (GlobalV::NSPIN == 2)
        {
            zipname = "output_DM1.npz";
            this->read_mat_npz(zipname, *(dm->get_DMR_pointer(2)));
        }

        this->pelec->psiToRho(*this->psi);

        int nspin0 = GlobalV::NSPIN == 2 ? 2 : 1;
        for (int is = 0; is < nspin0; is++)
        {
            std::string fn = GlobalV::global_out_dir + "/SPIN" + std::to_string(is + 1) + "_CHG.cube";
            ModuleIO::write_cube(
#ifdef __MPI
                this->pw_big->bz,
                this->pw_big->nbz,
                this->pw_rhod->nplane,
                this->pw_rhod->startz_current,
#endif
                this->pelec->charge->rho[is],
                is,
                GlobalV::NSPIN,
                istep,
                fn,
                this->pw_rhod->nx,
                this->pw_rhod->ny,
                this->pw_rhod->nz,
                this->pelec->eferm.get_efval(is),
                &(GlobalC::ucell),
                3,
                1);
        }

        return;
    }

    // the electron charge density should be symmetrized,
    // here is the initialization
    Symmetry_rho srho;
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        srho.begin(is, *(this->pelec->charge), this->pw_rho, GlobalC::Pgrid, GlobalC::ucell.symm);
    }

    // 1. calculate ewald energy.
    // mohan update 2021-02-25
    if (!PARAM.inp.test_skip_ewald)
    {
        this->pelec->f_en.ewald_energy = H_Ewald_pw::compute_ewald(GlobalC::ucell, this->pw_rho, this->sf.strucFac);
    }

    this->p_hamilt->non_first_scf = istep;
    return;
}


template class ESolver_KS_LCAO<double, double>;
template class ESolver_KS_LCAO<std::complex<double>, double>;
template class ESolver_KS_LCAO<std::complex<double>, std::complex<double>>;
} // namespace ModuleESolver
