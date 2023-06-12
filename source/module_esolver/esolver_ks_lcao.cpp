#include "esolver_ks_lcao.h"

#include "module_io/dos_nao.h"
#include "module_io/mulliken_charge.h"
#include "module_io/nscf_band.h"
#include "module_io/write_HS.h"
#include "module_io/write_istate_info.h"
#include "module_io/write_proj_band_lcao.h"
#include "module_io/output_log.h"

//--------------temporary----------------------------
#include "module_base/global_function.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_elecstate/occupy.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/print_info.h"
#ifdef __EXX
// #include "module_rpa/rpa.h"
#include "module_ri/RPA_LRI.h"
#endif

#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#include "module_hamilt_lcao/module_deepks/LCAO_deepks_interface.h"
#endif
//-----force& stress-------------------
#include "module_hamilt_lcao/hamilt_lcaodft/FORCE_STRESS.h"

//-----HSolver ElecState Hamilt--------
#include "module_elecstate/elecstate_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hsolver/hsolver_lcao.h"
// function used by deepks
#include "module_elecstate/cal_dm.h"
//---------------------------------------------------

namespace ModuleESolver
{

ESolver_KS_LCAO::ESolver_KS_LCAO()
{
    classname = "ESolver_KS_LCAO";
    basisname = "LCAO";
#ifdef __EXX
    if (GlobalC::exx_info.info_ri.real_number)
    {
        this->exx_lri_double = std::make_shared<Exx_LRI<double>>(GlobalC::exx_info.info_ri);
        this->exd = std::make_shared<Exx_LRI_Interface<double>>(this->exx_lri_double);
        this->LM.Hexxd = &this->exd->get_Hexxs();
    }
    else
    {
        this->exx_lri_complex = std::make_shared<Exx_LRI<std::complex<double>>>(GlobalC::exx_info.info_ri);
        this->exc = std::make_shared<Exx_LRI_Interface<std::complex<double>>>(this->exx_lri_complex);
        this->LM.Hexxc = &this->exc->get_Hexxs();
    }
#endif
}
ESolver_KS_LCAO::~ESolver_KS_LCAO()
{
    this->orb_con.clear_after_ions(GlobalC::UOT, GlobalC::ORB, GlobalV::deepks_setorb, GlobalC::ucell.infoNL.nproj);
}

void ESolver_KS_LCAO::Init(Input& inp, UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "Init");
    // if we are only calculating S, then there is no need
    // to prepare for potentials and so on

    if (GlobalV::CALCULATION == "get_S")
    {
        ucell.read_pseudo(GlobalV::ofs_running);

        if (ModuleSymmetry::Symmetry::symm_flag == 1)
        {
            this->symm.analy_sys(ucell, GlobalV::ofs_running);
            ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SYMMETRY");
        }

        // Setup the k points according to symmetry.
        kv.set(this->symm, GlobalV::global_kpoint_card, GlobalV::NSPIN, ucell.G, ucell.latvec);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT K-POINTS");

        Print_Info::setup_parameters(ucell, kv);
    }
    else
    {
        ESolver_KS::Init(inp, ucell);
    } // end ifnot get_S

    // init ElecState
    // autoset nbands in ElecState, it should before basis_init (for Psi 2d divid)
    if (this->pelec == nullptr)
    {
        this->pelec = new elecstate::ElecStateLCAO(&(chr),
                                                   &(kv),
                                                   kv.nks,
                                                   &(this->LOC),
                                                   &(this->UHM),
                                                   &(this->LOWF),
                                                   this->pw_rho,
                                                   pw_big);
    }

    //------------------init Basis_lcao----------------------
    // Init Basis should be put outside of Ensolver.
    // * reading the localized orbitals/projectors
    // * construct the interpolation tables.
    this->Init_Basis_lcao(this->orb_con, inp, ucell);
    //------------------init Basis_lcao----------------------

    // pass Hamilt-pointer to Operator
    this->UHM.genH.LM = this->UHM.LM = &this->LM;
    // pass basis-pointer to EState and Psi
    this->LOC.ParaV = this->LOWF.ParaV = this->LM.ParaV = &(this->orb_con.ParaV);

    if (GlobalV::CALCULATION == "get_S")
    {
        return;
    }

    //------------------init Hamilt_lcao----------------------
    // * allocate H and S matrices according to computational resources
    // * set the 'trace' between local H/S and global H/S
    this->LM.divide_HS_in_frag(GlobalV::GAMMA_ONLY_LOCAL, orb_con.ParaV, kv.nks);
    //------------------init Hamilt_lcao----------------------

#ifdef __EXX
    // PLEASE simplify the Exx_Global interface
    // mohan add 2021-03-25
    // Peize Lin 2016-12-03
    if (GlobalV::CALCULATION == "scf" || GlobalV::CALCULATION == "relax" || GlobalV::CALCULATION == "cell-relax"
        || GlobalV::CALCULATION == "md")
    {
        if (GlobalC::exx_info.info_global.cal_exx)
        {
            /* In the special "two-level" calculation case,
            first scf iteration only calculate the functional without exact exchange.
            but in "nscf" calculation, there is no need of "two-level" method. */
            if (ucell.atoms[0].ncpp.xc_func == "HSE" || ucell.atoms[0].ncpp.xc_func == "PBE0")
            {
                XC_Functional::set_xc_type("pbe");
            }
            else if (ucell.atoms[0].ncpp.xc_func == "SCAN0")
            {
                XC_Functional::set_xc_type("scan");
            }

            // GlobalC::exx_lcao.init();
            if (GlobalC::exx_info.info_ri.real_number)
                this->exx_lri_double->init(MPI_COMM_WORLD, kv);
            else
                this->exx_lri_complex->init(MPI_COMM_WORLD, kv);
        }
    }
#endif

    // Quxin added for DFT+U
    if (GlobalV::dft_plus_u)
    {
        GlobalC::dftu.init(ucell, this->LM, kv.nks);
    }

    // output is GlobalC::ppcell.vloc 3D local pseudopotentials
    // without structure factors
    // this function belongs to cell LOOP
    GlobalC::ppcell.init_vloc(GlobalC::ppcell.vloc, pw_rho);

    // init HSolver
    if (this->phsol == nullptr)
    {
        this->phsol = new hsolver::HSolverLCAO(this->LOWF.ParaV);
        this->phsol->method = GlobalV::KS_SOLVER;
    }

    // Inititlize the charge density.
    this->pelec->charge->allocate(GlobalV::NSPIN);
    this->pelec->omega = GlobalC::ucell.omega;

    // Initialize the potential.
    if (this->pelec->pot == nullptr)
    {
        this->pelec->pot = new elecstate::Potential(pw_rho,
                                                    &GlobalC::ucell,
                                                    &(GlobalC::ppcell.vloc),
                                                    &(sf),
                                                    &(this->pelec->f_en.etxc),
                                                    &(this->pelec->f_en.vtxc));
    }

#ifdef __DEEPKS
    // wenfei 2021-12-19
    // if we are performing DeePKS calculations, we need to load a model
    if (GlobalV::deepks_scf)
    {
        // load the DeePKS model from deep neural network
        GlobalC::ld.load_model(INPUT.deepks_model);
    }
#endif

    // Fix pelec->wg by ocp_kb
    if (GlobalV::ocp)
    {
        this->pelec->fixed_weights(GlobalV::ocp_kb);
    }
}

void ESolver_KS_LCAO::init_after_vc(Input& inp, UnitCell& ucell)
{
    ESolver_KS::init_after_vc(inp, ucell);

    if (GlobalV::md_prec_level == 2)
    {
        delete this->pelec;
        this->pelec = new elecstate::ElecStateLCAO(&(chr),
                                                   &(kv),
                                                   kv.nks,
                                                   &(this->LOC),
                                                   &(this->UHM),
                                                   &(this->LOWF),
                                                   this->pw_rho,
                                                   pw_big);

        GlobalC::ppcell.init_vloc(GlobalC::ppcell.vloc, pw_rho);

        this->pelec->charge->allocate(GlobalV::NSPIN);
        this->pelec->omega = GlobalC::ucell.omega;

        if (this->pelec->pot != nullptr)
        {
            delete this->pelec->pot;
            this->pelec->pot = new elecstate::Potential(pw_rho,
                                                        &GlobalC::ucell,
                                                        &(GlobalC::ppcell.vloc),
                                                        &(sf),
                                                        &(this->pelec->f_en.etxc),
                                                        &(this->pelec->f_en.vtxc));
        }
    }
}

double ESolver_KS_LCAO::cal_Energy()
{
    return this->pelec->f_en.etot;
}

void ESolver_KS_LCAO::cal_Force(ModuleBase::matrix& force)
{
    Force_Stress_LCAO FSL(this->RA, GlobalC::ucell.nat);
    FSL.getForceStress(GlobalV::CAL_FORCE,
                       GlobalV::CAL_STRESS,
                       GlobalV::TEST_FORCE,
                       GlobalV::TEST_STRESS,
                       this->LOC,
                       this->pelec,
                       this->psid,
                       this->psi,
                       this->UHM,
                       force,
                       this->scs,
                       this->sf,
                       this->kv,
                       this->pw_rho,
#ifdef __EXX
                        *this->exx_lri_double,
                        *this->exx_lri_complex,
#endif  
                       &this->symm);
    // delete RA after cal_Force
    this->RA.delete_grid();
    this->have_force = true;
}

void ESolver_KS_LCAO::cal_Stress(ModuleBase::matrix& stress)
{
    if (!this->have_force)
    {
        ModuleBase::matrix fcs;
        this->cal_Force(fcs);
    }
    stress = this->scs; // copy the stress
    this->have_force = false;
}

void ESolver_KS_LCAO::postprocess()
{
    GlobalV::ofs_running << "\n\n --------------------------------------------" << std::endl;
    GlobalV::ofs_running << std::setprecision(16);
    GlobalV::ofs_running << " !FINAL_ETOT_IS " << this->pelec->f_en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl;
    GlobalV::ofs_running << " --------------------------------------------\n\n" << std::endl;

    if (INPUT.out_dos != 0 || INPUT.out_band != 0 || INPUT.out_proj_band != 0)
    {
        GlobalV::ofs_running << "\n\n\n\n";
        GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
        GlobalV::ofs_running << " |                                                                    |" << std::endl;
        GlobalV::ofs_running << " | Post-processing of data:                                           |" << std::endl;
        GlobalV::ofs_running << " | DOS (density of states) and bands will be output here.             |" << std::endl;
        GlobalV::ofs_running << " | If atomic orbitals are used, Mulliken charge analysis can be done. |" << std::endl;
        GlobalV::ofs_running << " | Also the .bxsf file containing fermi surface information can be    |" << std::endl;
        GlobalV::ofs_running << " | done here.                                                         |" << std::endl;
        GlobalV::ofs_running << " |                                                                    |" << std::endl;
        GlobalV::ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
        GlobalV::ofs_running << "\n\n\n\n";
    }
    // qianrui modify 2020-10-18
    if (GlobalV::CALCULATION == "scf" || GlobalV::CALCULATION == "md" || GlobalV::CALCULATION == "relax")
    {
        ModuleIO::write_istate_info(this->pelec->ekb, this->pelec->wg, kv, &(GlobalC::Pkpoints));
    }

    int nspin0 = (GlobalV::NSPIN == 2) ? 2 : 1;

    if (INPUT.out_band) // pengfei 2014-10-13
    {
        int nks = 0;
        if (nspin0 == 1)
        {
            nks = kv.nkstot;
        }
        else if (nspin0 == 2)
        {
            nks = kv.nkstot / 2;
        }

        for (int is = 0; is < nspin0; is++)
        {
            std::stringstream ss2;
            ss2 << GlobalV::global_out_dir << "BANDS_" << is + 1 << ".dat";
            GlobalV::ofs_running << "\n Output bands in file: " << ss2.str() << std::endl;
            ModuleIO::nscf_band(is, ss2.str(), nks, GlobalV::NBANDS, 0.0, this->pelec->ekb, kv, &(GlobalC::Pkpoints));
        }
    } // out_band

    if (INPUT.out_proj_band) // Projeced band structure added by jiyy-2022-4-20
    {
        ModuleIO::write_proj_band_lcao(this->psid,
                                       this->psi,
                                       this->UHM,
                                       this->pelec,
                                       kv,
                                       GlobalC::ucell,
                                       GlobalC::ORB,
                                       GlobalC::GridD);
    }

    if (INPUT.out_dos)
    {
        ModuleIO::out_dos_nao(this->psid,
                               this->psi,
                               this->UHM,
                               this->pelec->ekb,
                               this->pelec->wg,
                               INPUT.dos_edelta_ev,
                               INPUT.dos_scale,
                               INPUT.dos_sigma,
                               *(this->pelec->klist),
                               GlobalC::Pkpoints,
                               GlobalC::ucell,
                               this->pelec->eferm,
                               GlobalV::NBANDS);
    }
}

void ESolver_KS_LCAO::Init_Basis_lcao(ORB_control& orb_con, Input& inp, UnitCell& ucell)
{
    // Set the variables first
    this->orb_con.gamma_only = GlobalV::GAMMA_ONLY_LOCAL;
    this->orb_con.nlocal = GlobalV::NLOCAL;
    this->orb_con.nbands = GlobalV::NBANDS;
    this->orb_con.ParaV.nspin = GlobalV::NSPIN;
    this->orb_con.dsize = GlobalV::DSIZE;
    this->orb_con.nb2d = GlobalV::NB2D;
    this->orb_con.dcolor = GlobalV::DCOLOR;
    this->orb_con.drank = GlobalV::DRANK;
    this->orb_con.myrank = GlobalV::MY_RANK;
    this->orb_con.calculation = GlobalV::CALCULATION;
    this->orb_con.ks_solver = GlobalV::KS_SOLVER;
    this->orb_con.setup_2d = true;

    // * reading the localized orbitals/projectors
    // * construct the interpolation tables.
    this->orb_con.read_orb_first(GlobalV::ofs_running,
                                 GlobalC::ORB,
                                 ucell.ntype,
                                 GlobalV::global_orbital_dir,
                                 ucell.orbital_fn,
                                 ucell.descriptor_file,
                                 ucell.lmax,
                                 inp.lcao_ecut,
                                 inp.lcao_dk,
                                 inp.lcao_dr,
                                 inp.lcao_rmax,
                                 GlobalV::deepks_setorb,
                                 inp.out_mat_r,
                                 GlobalV::CAL_FORCE,
                                 GlobalV::MY_RANK);

    ucell.infoNL.setupNonlocal(ucell.ntype, ucell.atoms, GlobalV::ofs_running, GlobalC::ORB);

    int Lmax = 0;
#ifdef __EXX
    Lmax = GlobalC::exx_info.info_ri.abfs_Lmax;
#endif
    this->orb_con.set_orb_tables(GlobalV::ofs_running,
                                 GlobalC::UOT,
                                 GlobalC::ORB,
                                 ucell.lat0,
                                 GlobalV::deepks_setorb,
                                 Lmax,
                                 ucell.infoNL.nprojmax,
                                 ucell.infoNL.nproj,
                                 ucell.infoNL.Beta);

    if (this->orb_con.setup_2d)
        this->orb_con.setup_2d_division(GlobalV::ofs_running, GlobalV::ofs_warning);
}

void ESolver_KS_LCAO::eachiterinit(const int istep, const int iter)
{

    // mohan add 2010-07-16
    // used for pulay mixing.
    if (iter == 1)
        this->p_chgmix->reset();

    // mohan update 2012-06-05
    this->pelec->f_en.deband_harris = this->pelec->cal_delta_eband();

    // mohan move it outside 2011-01-13
    // first need to calculate the weight according to
    // electrons number.

    if (wf.init_wfc == "file")
    {
        if (iter == 1)
        {
            std::cout << " WAVEFUN -> CHARGE " << std::endl;

            // calculate the density matrix using read in wave functions
            // and the ncalculate the charge density on grid.

            if (this->psi != nullptr)
            {
                this->pelec->psiToRho(this->psi[0]);
            }
            else
            {
                this->pelec->psiToRho(this->psid[0]);
            }

            // calculate the local potential(rho) again.
            // the grid integration will do in later grid integration.

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // a puzzle remains here.
            // if I don't renew potential,
            // The scf_thr is very small.
            // OneElectron, Hartree and
            // Exc energy are all correct
            // except the band energy.
            //
            // solved by mohan 2010-09-10
            // there are there rho here:
            // rho1: formed by read in orbitals.
            // rho2: atomic rho, used to construct H
            // rho3: generated by after diagonalize
            // here converged because rho3 and rho1
            // are very close.
            // so be careful here, make sure
            // rho1 and rho2 are the same rho.
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if (GlobalV::NSPIN == 4)
                GlobalC::ucell.cal_ux();
            this->pelec->pot->update_from_charge(this->pelec->charge, &GlobalC::ucell);
            this->pelec->f_en.descf = this->pelec->cal_delta_escf();
        }
    }

#ifdef __EXX
    // calculate exact-exchange
    if (GlobalC::exx_info.info_ri.real_number)
        this->exd->exx_eachiterinit(this->LOC, *(this->p_chgmix), iter);
    else
        this->exc->exx_eachiterinit(this->LOC, *(this->p_chgmix), iter);
#endif

    if (GlobalV::dft_plus_u)
    {
        GlobalC::dftu.cal_slater_UJ(pelec->charge->rho, pw_rho->nrxx); // Calculate U and J if Yukawa potential is used
    }

#ifdef __DEEPKS
    // the density matrixes of DeePKS have been updated in each iter
    GlobalC::ld.set_hr_cal(true);
#endif

    if (!GlobalV::GAMMA_ONLY_LOCAL)
    {
        this->UHM.GK.renew();
    }
}
void ESolver_KS_LCAO::hamilt2density(int istep, int iter, double ethr)
{
    // save input rho
    pelec->charge->save_rho_before_sum_band();

    // using HSolverLCAO::solve()
    if (this->phsol != nullptr)
    {
        // reset energy
        this->pelec->f_en.eband = 0.0;
        this->pelec->f_en.demet = 0.0;
        if (this->psi != nullptr)
        {
            this->phsol->solve(this->p_hamilt, this->psi[0], this->pelec, GlobalV::KS_SOLVER);
        }
        else if (this->psid != nullptr)
        {
            this->phsol->solve(this->p_hamilt, this->psid[0], this->pelec, GlobalV::KS_SOLVER);
        }

        if (GlobalV::out_bandgap)
        {
            if (!GlobalV::TWO_EFERMI)
            {
                this->pelec->cal_bandgap();
            }
            else
            {
                this->pelec->cal_bandgap_updw();
            }
        }
    }
    else
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW", "HSolver has not been initialed!");
    }

    // print ekb for each k point and each band
    for (int ik = 0; ik < kv.nks; ++ik)
    {
        this->pelec->print_band(ik, INPUT.printe, iter);
    }

#ifdef __EXX
    if (GlobalC::exx_info.info_ri.real_number)
        this->exd->exx_hamilt2density(*this->pelec, *this->LOWF.ParaV);
    else
        this->exc->exx_hamilt2density(*this->pelec, *this->LOWF.ParaV);
#endif

    // if DFT+U calculation is needed, this function will calculate
    // the local occupation number matrix and energy correction
    if (GlobalV::dft_plus_u)
    {
        if (GlobalC::dftu.omc != 2)
        {
            if (GlobalV::GAMMA_ONLY_LOCAL)
                GlobalC::dftu.cal_occup_m_gamma(iter, this->LOC.dm_gamma, this->p_chgmix->get_mixing_beta());
            else
                GlobalC::dftu.cal_occup_m_k(iter, this->LOC.dm_k, kv, this->p_chgmix->get_mixing_beta());
        }
        GlobalC::dftu.cal_energy_correction(istep);
        GlobalC::dftu.output();
    }

#ifdef __DEEPKS
    if (GlobalV::deepks_scf)
    {
        const Parallel_Orbitals* pv = this->LOWF.ParaV;
        if (GlobalV::GAMMA_ONLY_LOCAL)
        {
            GlobalC::ld.cal_e_delta_band(this->LOC.dm_gamma, pv->trace_loc_row, pv->trace_loc_col, pv->nrow);
        }
        else
        {
            GlobalC::ld
                .cal_e_delta_band_k(this->LOC.dm_k, pv->trace_loc_row, pv->trace_loc_col, kv.nks, pv->nrow, pv->ncol);
        }
    }
#endif
    // (4) mohan add 2010-06-24
    // using new charge density.
    this->pelec->cal_energies(1);

    // (5) symmetrize the charge density
    Symmetry_rho srho;
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        srho.begin(is, *(pelec->charge), pw_rho, GlobalC::Pgrid, this->symm);
    }

    // (6) compute magnetization, only for spin==2
    GlobalC::ucell.magnet.compute_magnetization(this->pelec->charge->nrxx,
                                                this->pelec->charge->nxyz,
                                                this->pelec->charge->rho,
                                                this->pelec->nelec_spin.data());

    // (7) calculate delta energy
    this->pelec->f_en.deband = this->pelec->cal_delta_eband();
}
void ESolver_KS_LCAO::updatepot(const int istep, const int iter)
{
    // print Hamiltonian and Overlap matrix
    if (this->conv_elec)
    {
        if (!GlobalV::GAMMA_ONLY_LOCAL && hsolver::HSolverLCAO::out_mat_hs)
        {
            this->UHM.GK.renew(true);
        }
        for (int ik = 0; ik < kv.nks; ++ik)
        {
            if (hsolver::HSolverLCAO::out_mat_hs)
            {
                this->p_hamilt->updateHk(ik);
            }
            bool bit = false; // LiuXh, 2017-03-21
            // if set bit = true, there would be error in soc-multi-core calculation, noted by zhengdy-soc
            if (this->psi != nullptr)
            {
                hamilt::MatrixBlock<complex<double>> h_mat, s_mat;
                this->p_hamilt->matrix(h_mat, s_mat);
                ModuleIO::saving_HS(h_mat.p,
                                    s_mat.p,
                                    bit,
                                    hsolver::HSolverLCAO::out_mat_hs,
                                    "data-" + std::to_string(ik),
                                    this->LOWF.ParaV[0],
                                    1); // LiuXh, 2017-03-21
            }
            else if (this->psid != nullptr)
            { // gamma_only case, Hloc and Sloc are correct H and S matrix
                ModuleIO::saving_HS(this->LM.Hloc.data(),
                                    this->LM.Sloc.data(),
                                    bit,
                                    hsolver::HSolverLCAO::out_mat_hs,
                                    "data-" + std::to_string(ik),
                                    this->LOWF.ParaV[0],
                                    1); // LiuXh, 2017-03-21
            }
        }
    }

    if (this->conv_elec)
    {
        if (elecstate::ElecStateLCAO::out_wfc_lcao)
        {
            elecstate::ElecStateLCAO::out_wfc_flag = 1;
        }
        for (int ik = 0; ik < kv.nks; ik++)
        {
            if (this->psi != nullptr)
            {
                this->psi[0].fix_k(ik);
                this->pelec->print_psi(this->psi[0]);
            }
            else
            {
                this->psid[0].fix_k(ik);
                this->pelec->print_psi(this->psid[0]);
            }
        }
        elecstate::ElecStateLCAO::out_wfc_flag = 0;
    }

    // (9) Calculate new potential according to new Charge Density.

    if (this->conv_elec || iter == GlobalV::SCF_NMAX)
    {
        if (GlobalV::out_pot < 0) // mohan add 2011-10-10
        {
            GlobalV::out_pot = -2;
        }
    }
    if (!this->conv_elec)
    {
        if (GlobalV::NSPIN == 4)
            GlobalC::ucell.cal_ux();
        this->pelec->pot->update_from_charge(this->pelec->charge, &GlobalC::ucell);
        this->pelec->f_en.descf = this->pelec->cal_delta_escf();
    }
    else
    {
        this->pelec->cal_converged();
    }
}
void ESolver_KS_LCAO::eachiterfinish(int iter)
{
    //-----------------------------------
    // save charge density
    //-----------------------------------
    // Peize Lin add 2020.04.04
    if (GlobalC::restart.info_save.save_charge)
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            GlobalC::restart.save_disk(*this->UHM.LM, "charge", is, pelec->charge->nrxx, pelec->charge->rho);
        }
    }

    //-----------------------------------
    // output charge density for tmp
    //-----------------------------------
    bool print = false;
    if (this->out_freq_elec && iter % this->out_freq_elec == 0)
    {
        print = true;
    }

    if (print)
    {
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            this->create_Output_Rho(is, iter, "tmp_").write();
            this->create_Output_DM(is, iter).write();
            if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
            {
                this->create_Output_Kin(is, iter, "tmp_").write();
            }
        }
    }

    // (11) calculate the total energy.
    this->pelec->cal_energies(2);
}

void ESolver_KS_LCAO::afterscf(const int istep)
{
    if (this->LOC.out_dm1 == 1)
    {
        this->create_Output_DM1(istep).write();
    }

    if (GlobalV::out_chg)
    {
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            this->create_Output_Rho(is, istep).write();
            if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
            {
                this->create_Output_Kin(is, istep).write();
            }
        }
    }

    if (this->LOC.out_dm)
    {
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            this->create_Output_DM(is, istep).write();
        }
    }

#ifdef __EXX
    if (GlobalC::exx_info.info_global.cal_exx) // Peize Lin add if 2022.11.14
    {
        const std::string file_name_exx = GlobalV::global_out_dir + "HexxR_" + std::to_string(GlobalV::MY_RANK);
        if (GlobalC::exx_info.info_ri.real_number)
            this->exd->write_Hexxs(file_name_exx);
        else
            this->exc->write_Hexxs(file_name_exx);
    }
#endif

    this->create_Output_Potential(istep).write();

    ModuleIO::output_convergence_after_scf(this->conv_elec, this->pelec->f_en.etot);
    ModuleIO::output_efermi(this->conv_elec, this->pelec->eferm.ef);

    if (GlobalV::OUT_LEVEL != "m")
    {
        this->pelec->print_eigenvalue(GlobalV::ofs_running);
    }

#ifdef __DEEPKS
    std::shared_ptr<LCAO_Deepks> ld_shared_ptr(&GlobalC::ld,[](LCAO_Deepks*){});
    LCAO_Deepks_Interface LDI = LCAO_Deepks_Interface(ld_shared_ptr);
    LDI.out_deepks_labels(this->pelec->f_en.etot,
                          this->pelec->klist->nks,
                          GlobalC::ucell.nat,
                          this->pelec->ekb,
                          this->pelec->klist->kvec_d,
                          GlobalC::ucell,
                          GlobalC::ORB,
                          GlobalC::GridD,
                          this->LOWF.ParaV,
                          *(this->psi),
                          *(this->psid),
                          this->LOC.dm_gamma,
                          this->LOC.dm_k);

#endif
    // 3. some outputs
#ifdef __EXX
    if (INPUT.rpa)
    {
        // ModuleRPA::DFT_RPA_interface rpa_interface(GlobalC::exx_info.info_global);
        // rpa_interface.rpa_exx_lcao().info.files_abfs = GlobalV::rpa_orbitals;
        // rpa_interface.out_for_RPA(*(this->LOWF.ParaV), *(this->psi), this->LOC, this->pelec);
        RPA_LRI<double> rpa_lri_double(GlobalC::exx_info.info_ri);
        rpa_lri_double.cal_postSCF_exx(this->LOC, MPI_COMM_WORLD, kv, *this->LOWF.ParaV);
        rpa_lri_double.init(MPI_COMM_WORLD, kv);
        rpa_lri_double.out_for_RPA(*(this->LOWF.ParaV), *(this->psi), this->pelec);
    }
#endif

    this->create_Output_Mat_Sparse(istep).write();

    // GlobalV::mulliken charge analysis
    if (GlobalV::out_mul)
    {
        if (GlobalV::CALCULATION != "md" || (istep % GlobalV::out_interval == 0))
        {
            ModuleIO::out_mulliken(istep, this->UHM, this->LOC, kv);
        }
    } // qifeng add 2019/9/10, jiyy modify 2023/2/27, liuyu move here 2023-04-18

    if (!GlobalV::CAL_FORCE && !GlobalV::CAL_STRESS)
    {
        RA.delete_grid();
    }
}

bool ESolver_KS_LCAO::do_after_converge(int& iter)
{
#ifdef __EXX
    if (GlobalC::exx_info.info_ri.real_number)
        return this->exd->exx_after_converge(*this->p_hamilt, this->LM, this->LOC, kv, iter);
    else
        return this->exc->exx_after_converge(*this->p_hamilt, this->LM, this->LOC, kv, iter);
#endif // __EXX
    return true;
}

ModuleIO::Output_DM ESolver_KS_LCAO::create_Output_DM(int is, int iter)
{
    int precision = 3;
    return ModuleIO::Output_DM(this->GridT,
                               is,
                               iter,
                               precision,
                               this->LOC.out_dm,
                               this->LOC.DM,
                               this->pelec->eferm.get_efval(is),
                               &(GlobalC::ucell),
                               GlobalV::global_out_dir,
                               GlobalV::GAMMA_ONLY_LOCAL);
}

ModuleIO::Output_DM1 ESolver_KS_LCAO::create_Output_DM1(int istep)
{
    return ModuleIO::Output_DM1(GlobalV::NSPIN, istep, this->LOC, this->RA, this->kv);
}

ModuleIO::Output_Mat_Sparse ESolver_KS_LCAO::create_Output_Mat_Sparse(int istep)
{
    bool is_md = (GlobalV::CALCULATION == "md") ? true : false;
    return ModuleIO::Output_Mat_Sparse(hsolver::HSolverLCAO::out_mat_hsR,
                                       hsolver::HSolverLCAO::out_mat_t,
                                       hsolver::HSolverLCAO::out_mat_dh,
                                       INPUT.out_mat_r,
                                       is_md,
                                       istep,
                                       GlobalV::out_interval,
                                       this->pelec->pot->get_effective_v(),
                                       *this->LOWF.ParaV,
                                       this->UHM,
                                       this->LM,
                                       this->kv);
}

} // namespace ModuleESolver