#include "esolver_ks_lcao_tddft.h"

#include "src_io/cal_r_overlap_R.h"

//--------------temporary----------------------------
#include "../module_base/blas_connector.h"
#include "../module_base/global_function.h"
#include "../module_base/scalapack_connector.h"
#include "../src_io/print_info.h"
#include "../src_pw/global.h"
#include "input_update.h"
#include "src_io/chi0_hilbert.h"
#include "src_lcao/ELEC_evolve.h"
#include "src_pw/occupy.h"
#include "src_pw/symmetry_rho.h"
#include "src_pw/threshold_elec.h"

//-----HSolver ElecState Hamilt--------
#include "module_elecstate/elecstate_lcao.h"
#include "module_elecstate/elecstate_lcao_tddft.h"
#include "module_hamilt/hamilt_lcao.h"
#include "module_hsolver/hsolver_lcao.h"
#include "module_psi/psi.h"

//-----force& stress-------------------
#include "src_lcao/FORCE_STRESS.h"
#include "src_lcao/FORCE_STRESS_tddft.h"

//---------------------------------------------------

namespace ModuleESolver
{

ESolver_KS_LCAO_TDDFT::ESolver_KS_LCAO_TDDFT()
{
    classname = "ESolver_KS_LCAO_TDDFT";
    basisname = "LCAO";
}
ESolver_KS_LCAO_TDDFT::~ESolver_KS_LCAO_TDDFT()
{
    // this->orb_con.clear_after_ions(GlobalC::UOT, GlobalC::ORB, GlobalV::deepks_setorb, GlobalC::ucell.infoNL.nproj);
    delete psi_laststep;
    delete pelec_td;
}

void ESolver_KS_LCAO_TDDFT::Init(Input& inp, UnitCell_pseudo& ucell)
{
    ESolver_KS::Init(inp, ucell);

    // Inititlize the charge density.
    GlobalC::CHR.allocate(GlobalV::NSPIN, GlobalC::rhopw->nrxx, GlobalC::rhopw->npw);
    // GlobalC::CHR.allocate(GlobalV::NSPIN, GlobalC::rhopw->nrxx, GlobalC::rhopw->npw);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT CHARGE");
    // Initializee the potential.
    GlobalC::pot.allocate(GlobalC::rhopw->nrxx);

    // Initialize the local wave functions.
    // npwx, eigenvalues, and weights
    // npwx may change according to cell change
    // this function belongs to cell LOOP
    GlobalC::wf.allocate_ekb_wg(GlobalC::kv.nks);

    // Initialize the FFT.
    // this function belongs to cell LOOP

    // output is GlobalC::ppcell.vloc 3D local pseudopotentials
    // without structure factors
    // this function belongs to cell LOOP
    GlobalC::ppcell.init_vloc(GlobalC::ppcell.vloc, GlobalC::rhopw);

    // Initialize the sum of all local potentials.
    // if ion_step==0, read in/initialize the potentials
    // this function belongs to ions LOOP
    int ion_step = 0;
    GlobalC::pot.init_pot(ion_step, GlobalC::sf.strucFac);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT POTENTIAL");

    //------------------init Basis_lcao----------------------
    // Init Basis should be put outside of Ensolver.
    // * reading the localized orbitals/projectors
    // * construct the interpolation tables.
    this->Init_Basis_lcao(this->orb_con, inp, ucell);
    //------------------init Basis_lcao----------------------

    //------------------init Hamilt_lcao----------------------
    // * allocate H and S matrices according to computational resources
    // * set the 'trace' between local H/S and global H/S
    this->LM.divide_HS_in_frag(GlobalV::GAMMA_ONLY_LOCAL, orb_con.ParaV);
    //------------------init Hamilt_lcao----------------------

    // pass Hamilt-pointer to Operator
    this->UHM.genH.LM = this->UHM.LM = &this->LM;
    // pass basis-pointer to EState and Psi
    this->LOC.ParaV = this->LOWF.ParaV = this->LM.ParaV;

    // init Psi, HSolver, ElecState, Hamilt
    if (this->phsol != nullptr)
    {
        if (this->phsol->classname != "HSolverLCAO")
        {
            delete this->phsol;
            this->phsol = nullptr;
        }
    }
    else
    {
        this->phsol = new hsolver::HSolverLCAO(this->LOWF.ParaV);
        this->phsol->method = GlobalV::KS_SOLVER;
    }
    if (this->pelec_td != nullptr)
    {
        if (this->pelec_td->classname != "ElecStateLCAO")
        {
            delete this->pelec_td;
            this->pelec_td = nullptr;
        }
    }
    else
    {
        this->pelec_td = new elecstate::ElecStateLCAO_TDDFT((Charge*)(&(GlobalC::CHR)),
                                                            &(GlobalC::kv),
                                                            GlobalC::kv.nks,
                                                            GlobalV::NBANDS,
                                                            &(this->LOC),
                                                            &(this->UHM),
                                                            &(this->LOWF));
    }
}

void ESolver_KS_LCAO_TDDFT::eachiterinit(const int istep, const int iter)
{
    std::string ufile = "CHANGE";
    Update_input UI;
    UI.init(ufile);

    // mohan add 2010-07-16
    // used for pulay mixing.
    if (iter == 1)
    {
        GlobalC::CHR.set_new_e_iteration(true);
    }
    else
    {
        GlobalC::CHR.set_new_e_iteration(false);
    }

    if (GlobalV::FINAL_SCF && iter == 1)
    {
        GlobalC::CHR.irstep = 0;
        GlobalC::CHR.idstep = 0;
        GlobalC::CHR.totstep = 0;
    }

    // mohan update 2012-06-05
    GlobalC::en.calculate_harris(1);

    // mohan move it outside 2011-01-13
    // first need to calculate the weight according to
    // electrons number.
    // mohan add iter > 1 on 2011-04-02
    // because the GlobalC::en.ekb has not value now.
    // so the smearing can not be done.
    // if (iter > 1 && istep <= 1 && GlobalV::ocp == 0)
    //    Occupy::calculate_weights();

    if (GlobalC::wf.init_wfc == "file")
    {
        if (iter == 1)
        {
            std::cout << " WAVEFUN -> CHARGE " << std::endl;

            // The occupation should be read in together.
            // Occupy::calculate_weights(); //mohan add 2012-02-15

            // calculate the density matrix using read in wave functions
            // and the ncalculate the charge density on grid.

            // transform wg and ekb to elecstate first
            for (int ik = 0; ik < this->pelec_td->ekb.nr; ++ik)
            {
                for (int ib = 0; ib < this->pelec_td->ekb.nc; ++ib)
                {
                    this->pelec_td->ekb(ik, ib) = GlobalC::wf.ekb[ik][ib];
                    this->pelec_td->wg(ik, ib) = GlobalC::wf.wg(ik, ib);
                }
            }
            if (this->psi != nullptr)
            {
                if (istep >= 2)
                {
                    this->pelec_td->psiToRho_td(this->psi[0]);
                }
                else
                {
                    this->pelec_td->psiToRho(this->psi[0]);
                }
            }
            else
            {
                this->pelec_td->psiToRho(this->psid[0]);
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
            GlobalC::pot.vr = GlobalC::pot.v_of_rho(GlobalC::CHR.rho, GlobalC::CHR.rho_core);
            GlobalC::en.delta_escf();
            if (ELEC_evolve::td_vext == 0)
            {
                GlobalC::pot.set_vr_eff();
            }
            else
            {
                GlobalC::pot.set_vrs_tddft(istep);
            }
        }
    }

    if (!GlobalV::GAMMA_ONLY_LOCAL)
    {
        if (this->UHM.GK.get_spin() != -1)
        {
            int start_spin = -1;
            this->UHM.GK.reset_spin(start_spin);
            this->UHM.GK.destroy_pvpR();
            this->UHM.GK.allocate_pvpR();
        }
    }
}

void ESolver_KS_LCAO_TDDFT::hamilt2density(int istep, int iter, double ethr, ModuleBase::Vector3<double>* vel)
{

    GlobalC::CHR.save_rho_before_sum_band();

    if (ELEC_evolve::tddft && istep >= 2 && !GlobalV::GAMMA_ONLY_LOCAL)
    {
        ELEC_evolve::evolve_psi(istep,
                                this->p_hamilt,
                                this->LOWF,
                                this->psi,
                                this->psi_laststep,
                                this->RA,
                                this->UHM,
                                vel);
        this->pelec_td->psiToRho_td(this->psi[0]);
        // this->pelec_td->psiToRho(this->psi[0]);
    }
    // using HSolverLCAO::solve()
    else if (this->phsol != nullptr)
    {
        // reset energy
        this->pelec_td->eband = 0.0;
        this->pelec_td->demet = 0.0;
        this->pelec_td->ef = 0.0;
        GlobalC::en.ef_up = 0.0;
        GlobalC::en.ef_dw = 0.0;
        if (this->psi != nullptr)
        {
            this->phsol->solve(this->p_hamilt, this->psi[0], this->pelec_td, GlobalV::KS_SOLVER);
        }
        else if (this->psid != nullptr)
        {
            this->phsol->solve(this->p_hamilt, this->psid[0], this->pelec_td, GlobalV::KS_SOLVER);
        }
    }
    else
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_LCAO", "HSolver has not been initialed!");
    }

    if (iter == 1)
    {
        GlobalV::ofs_running
            << "------------------------------------------------------------------------------------------------"
            << endl;
        GlobalV::ofs_running << "occupation : ";
        for (int ik = 0; ik < GlobalC::kv.nks; ik++)
        {
            for (int ib = 0; ib < GlobalV::NBANDS; ib++)
            {
                GlobalV::ofs_running << this->pelec_td->wg(ik, ib) << " ";
            }
        }
        GlobalV::ofs_running << endl;
        GlobalV::ofs_running
            << "------------------------------------------------------------------------------------------------"
            << endl;
    }

    // transform energy for print
    GlobalC::en.eband = this->pelec_td->eband;
    GlobalC::en.demet = this->pelec_td->demet;
    GlobalC::en.ef = this->pelec_td->ef;

    // (3) sum bands to calculate charge density
    // if (istep <= 1 ) Occupy::calculate_weights();

    for (int ik = 0; ik < GlobalC::kv.nks; ++ik)
    {
        GlobalC::en.print_band(ik);
    }

    // (4) mohan add 2010-06-24
    // using new charge density.
    GlobalC::en.calculate_harris(2);

    // (5) symmetrize the charge density
    if (istep <= 1)
    {
        Symmetry_rho srho;
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            srho.begin(is, GlobalC::CHR, GlobalC::rhopw, GlobalC::Pgrid, GlobalC::symm);
        }
    }

    // (6) compute magnetization, only for spin==2
    GlobalC::ucell.magnet.compute_magnetization();

    // resume codes!
    //-------------------------------------------------------------------------
    // this->GlobalC::LOWF.init_Cij( 0 ); // check the orthogonality of local orbital.
    // GlobalC::CHR.sum_band(); use local orbital in plane wave basis to calculate bands.
    // but must has evc first!
    //-------------------------------------------------------------------------

    // (7) calculate delta energy
    GlobalC::en.deband = GlobalC::en.delta_e();
}

void ESolver_KS_LCAO_TDDFT::hamilt2density(int istep, int iter, double ethr)
{

    GlobalC::CHR.save_rho_before_sum_band();

    if (this->phsol != nullptr)
    {
        // reset energy
        this->pelec_td->eband = 0.0;
        this->pelec_td->demet = 0.0;
        this->pelec_td->ef = 0.0;
        GlobalC::en.ef_up = 0.0;
        GlobalC::en.ef_dw = 0.0;
        if (this->psi != nullptr)
        {
            this->phsol->solve(this->p_hamilt, this->psi[0], this->pelec_td, GlobalV::KS_SOLVER);
        }
        else if (this->psid != nullptr)
        {
            this->phsol->solve(this->p_hamilt, this->psid[0], this->pelec_td, GlobalV::KS_SOLVER);
        }
    }
    else
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_LCAO", "HSolver has not been initialed!");
    }

    if (iter == 1)
    {
        GlobalV::ofs_running
            << "------------------------------------------------------------------------------------------------"
            << endl;
        GlobalV::ofs_running << "occupation : ";
        for (int ik = 0; ik < GlobalC::kv.nks; ik++)
        {
            for (int ib = 0; ib < GlobalV::NBANDS; ib++)
            {
                GlobalV::ofs_running << this->pelec_td->wg(ik, ib) << " ";
            }
        }
        GlobalV::ofs_running << endl;
        GlobalV::ofs_running
            << "------------------------------------------------------------------------------------------------"
            << endl;
    }

    // transform energy for print
    GlobalC::en.eband = this->pelec_td->eband;
    GlobalC::en.demet = this->pelec_td->demet;
    GlobalC::en.ef = this->pelec_td->ef;

    // (3) sum bands to calculate charge density
    // if (istep <= 1 ) Occupy::calculate_weights();

    for (int ik = 0; ik < GlobalC::kv.nks; ++ik)
    {
        GlobalC::en.print_band(ik);
    }

    // (4) mohan add 2010-06-24
    // using new charge density.
    GlobalC::en.calculate_harris(2);

    // (5) symmetrize the charge density
    if (istep <= 1)
    {
        Symmetry_rho srho;
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            srho.begin(is, GlobalC::CHR, GlobalC::rhopw, GlobalC::Pgrid, GlobalC::symm);
        }
    }

    // (6) compute magnetization, only for spin==2
    GlobalC::ucell.magnet.compute_magnetization();

    // resume codes!
    //-------------------------------------------------------------------------
    // this->GlobalC::LOWF.init_Cij( 0 ); // check the orthogonality of local orbital.
    // GlobalC::CHR.sum_band(); use local orbital in plane wave basis to calculate bands.
    // but must has evc first!
    //-------------------------------------------------------------------------

    // (7) calculate delta energy
    GlobalC::en.deband = GlobalC::en.delta_e();
}

void ESolver_KS_LCAO_TDDFT::updatepot(const int istep, const int iter)
{
    // (9) Calculate new potential according to new Charge Density.

    if (this->conv_elec || iter == GlobalV::SCF_NMAX)
    {
        if (GlobalC::pot.out_pot < 0) // mohan add 2011-10-10
        {
            GlobalC::pot.out_pot = -2;
        }
    }
    if (!this->conv_elec)
    {
        GlobalC::pot.vr = GlobalC::pot.v_of_rho(GlobalC::CHR.rho, GlobalC::CHR.rho_core);
        GlobalC::en.delta_escf();
    }
    else
    {
        GlobalC::pot.vnew = GlobalC::pot.v_of_rho(GlobalC::CHR.rho, GlobalC::CHR.rho_core);
        //(used later for scf correction to the forces )
        GlobalC::pot.vnew -= GlobalC::pot.vr;
        GlobalC::en.descf = 0.0;
    }

    // add Vloc to Vhxc.
    if (ELEC_evolve::td_vext == 0)
    {
        GlobalC::pot.set_vr_eff();
    }
    else
    {
        GlobalC::pot.set_vrs_tddft(istep);
    }

    // store wfc
    if (istep >= 1 && this->conv_elec)
    {
        if (this->psi_laststep == nullptr)
#ifdef __MPI
            this->psi_laststep = new psi::Psi<std::complex<double>>(GlobalC::kv.nks,
                                                                    this->LOWF.ParaV->ncol_bands,
                                                                    this->LOWF.ParaV->nrow,
                                                                    nullptr);
#else
            this->psi_laststep
                = new psi::Psi<std::complex<double>>(GlobalC::kv.nks, GlobalV::NBANDS, GlobalV::NLOCAL, nullptr);
#endif
        std::complex<double>* tmp = psi[0].get_pointer();
        for (int index = 0; index < psi[0].size(); ++index)
            psi_laststep[0].get_pointer()[index] = tmp[index];
        if (istep > 1)
            this->cal_edm_tddft();
    }

    if (this->conv_elec)
    {
        GlobalV::ofs_running
            << "------------------------------------------------------------------------------------------------"
            << endl;
        GlobalV::ofs_running << "Eii : ";
        for (int ik = 0; ik < GlobalC::kv.nks; ik++)
        {
            for (int ib = 0; ib < GlobalV::NBANDS; ib++)
            {
                GlobalV::ofs_running << this->pelec_td->ekb(ik, ib) << " ";
            }
        }
        GlobalV::ofs_running << endl;
        GlobalV::ofs_running
            << "------------------------------------------------------------------------------------------------"
            << endl;
    }
}

void ESolver_KS_LCAO_TDDFT::afterscf(const int istep)
{
    for (int ik = 0; ik < this->pelec_td->ekb.nr; ++ik)
    {
        for (int ib = 0; ib < this->pelec_td->ekb.nc; ++ib)
        {
            GlobalC::wf.ekb[ik][ib] = this->pelec_td->ekb(ik, ib);
            GlobalC::wf.wg(ik, ib) = this->pelec_td->wg(ik, ib);
        }
    }
    // if (this->conv_elec || iter == GlobalV::SCF_NMAX)
    // {
    //--------------------------------------
    // 1. output charge density for converged,
    // 0 means don't need to consider iter,
    //--------------------------------------
    if (GlobalC::chi0_hilbert.epsilon) // pengfei 2016-11-23
    {
        std::cout << "eta = " << GlobalC::chi0_hilbert.eta << std::endl;
        std::cout << "domega = " << GlobalC::chi0_hilbert.domega << std::endl;
        std::cout << "nomega = " << GlobalC::chi0_hilbert.nomega << std::endl;
        std::cout << "dim = " << GlobalC::chi0_hilbert.dim << std::endl;
        // std::cout <<"oband = "<<GlobalC::chi0_hilbert.oband<<std::endl;
        GlobalC::chi0_hilbert.wfc_k_grid = this->LOWF.wfc_k_grid;
        GlobalC::chi0_hilbert.Chi();
    }

    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        const int precision = 3;

        std::stringstream ssc;
        std::stringstream ss1;
        ssc << GlobalV::global_out_dir << "SPIN" << is + 1 << "_CHG";
        ss1 << GlobalV::global_out_dir << "SPIN" << is + 1 << "_CHG.cube";
        GlobalC::CHR.write_rho(GlobalC::CHR.rho_save[is], is, 0, ssc.str()); // mohan add 2007-10-17
        GlobalC::CHR.write_rho_cube(GlobalC::CHR.rho_save[is], is, ss1.str(), 3);

        std::stringstream ssd;
        if (GlobalV::GAMMA_ONLY_LOCAL)
        {
            ssd << GlobalV::global_out_dir << "SPIN" << is + 1 << "_DM";
        }
        else
        {
            ssd << GlobalV::global_out_dir << "SPIN" << is + 1 << "_DM_R";
        }
        this->LOC.write_dm(is, 0, ssd.str(), precision);

        if (GlobalC::pot.out_pot == 1) // LiuXh add 20200701
        {
            std::stringstream ssp;
            ssp << GlobalV::global_out_dir << "SPIN" << is + 1 << "_POT";
            GlobalC::pot.write_potential(is, 0, ssp.str(), GlobalC::pot.vr_eff, precision);
        }
    }

    if (this->conv_elec)
    {
        GlobalV::ofs_running << "\n charge density convergence is achieved" << std::endl;
        GlobalV::ofs_running << " final etot is " << GlobalC::en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl;
    }

    if (GlobalV::OUT_LEVEL != "m")
    {
        // Threshold_Elec::print_eigenvalue(GlobalV::ofs_running);
    }

    if (this->conv_elec)
    {
        // xiaohui add "OUT_LEVEL", 2015-09-16
        if (GlobalV::OUT_LEVEL != "m")
            GlobalV::ofs_running << std::setprecision(16);
        if (GlobalV::OUT_LEVEL != "m")
            GlobalV::ofs_running << " EFERMI = " << GlobalC::en.ef * ModuleBase::Ry_to_eV << " eV" << std::endl;
        if (GlobalV::OUT_LEVEL == "ie")
        {
            GlobalV::ofs_running << " " << GlobalV::global_out_dir << " final etot is "
                                 << GlobalC::en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl;
        }
    }
    else
    {
        GlobalV::ofs_running << " !! convergence has not been achieved @_@" << std::endl;
        if (GlobalV::OUT_LEVEL == "ie" || GlobalV::OUT_LEVEL == "m") // xiaohui add "m" option, 2015-09-16
            std::cout << " !! CONVERGENCE HAS NOT BEEN ACHIEVED !!" << std::endl;
    }

    if (hsolver::HSolverLCAO::out_mat_hsR)
    {
        this->output_HS_R(istep); // LiuXh add 2019-07-15
    }

    // add by jingan for out r_R matrix 2019.8.14
    if (INPUT.out_mat_r)
    {
        cal_r_overlap_R r_matrix;
        r_matrix.init(*this->LOWF.ParaV);

        if (hsolver::HSolverLCAO::out_mat_hsR)
        {
            r_matrix.out_rR_other(this->LM.output_R_coor);
        }
        else
        {
            r_matrix.out_rR();
        }
    }
}

void ESolver_KS_LCAO_TDDFT::Run(const int istep, UnitCell_pseudo& ucell, ModuleBase::Vector3<double>* vel)
{
    if (!(GlobalV::CALCULATION == "scf" || GlobalV::CALCULATION == "md" || GlobalV::CALCULATION == "relax"
          || GlobalV::CALCULATION == "cell-relax" || GlobalV::CALCULATION.substr(0, 3) == "sto"))
    {
        this->othercalculation(istep);
    }
    else
    {
        ModuleBase::timer::tick(this->classname, "Run");

        if (this->maxniter > 0)
            this->printhead(); // print the headline on the screen.
        this->beforescf(istep); // Something else to do before the iter loop

        bool firstscf = true;
        this->conv_elec = false;
        this->niter = this->maxniter;
        for (int iter = 1; iter <= this->maxniter; ++iter)
        {
            writehead(GlobalV::ofs_running, istep, iter);
#ifdef __MPI
            auto iterstart = MPI_Wtime();
#else
            auto iterstart = std::chrono::system_clock::now();
#endif
            double diag_ethr = this->phsol->set_diagethr(istep, iter, drho);
            eachiterinit(istep, iter);
            this->hamilt2density(istep, iter, diag_ethr, vel);

            //<Temporary> It may be changed when more clever parallel algorithm is put forward.
            // When parallel algorithm for bands are adopted. Density will only be treated in the first group.
            //(Different ranks should have abtained the same, but small differences always exist in practice.)
            // Maybe in the future, density and wavefunctions should use different parallel algorithms, in which
            // they do not occupy all processors, for example wavefunctions uses 20 processors while density uses 10.
            if (GlobalV::MY_STOGROUP == 0)
            {
                // double drho = this->estate.caldr2();
                // EState should be used after it is constructed.
                drho = GlobalC::CHR.get_drho();
                double hsolver_error = 0.0;
                if (firstscf)
                {
                    firstscf = false;
                    hsolver_error = this->phsol->cal_hsolerror();
                    // The error of HSolver is larger than drho, so a more precise HSolver should be excuconv_elected.
                    if (hsolver_error > drho)
                    {
                        diag_ethr = this->phsol->reset_diagethr(GlobalV::ofs_running, hsolver_error, drho);
                        this->hamilt2density(istep, iter, diag_ethr, vel);
                        drho = GlobalC::CHR.get_drho();
                        hsolver_error = this->phsol->cal_hsolerror();
                    }
                }

                this->conv_elec = (drho < this->scf_thr);

                // If drho < hsolver_error in the first iter or drho < scf_thr, we do not change rho.
                if (drho < hsolver_error || this->conv_elec)
                {
                    if (drho < hsolver_error)
                        GlobalV::ofs_warning << " drho < hsolver_error, keep charge density unchanged." << std::endl;
                }
                else
                {
                    // charge mixing
                    // conv_elec = this->estate.mix_rho();
                    GlobalC::CHR.mix_rho(iter);
                }
            }
#ifdef __MPI
            MPI_Bcast(&drho, 1, MPI_DOUBLE, 0, PARAPW_WORLD);
            MPI_Bcast(&this->conv_elec, 1, MPI_DOUBLE, 0, PARAPW_WORLD);
            MPI_Bcast(GlobalC::CHR.rho[0], GlobalC::rhopw->nrxx, MPI_DOUBLE, 0, PARAPW_WORLD);
#endif

            // Hamilt should be used after it is constructed.
            // this->phamilt->update(conv_elec);
            updatepot(istep, iter);
            eachiterfinish(iter);
#ifdef __MPI
            double duration = (double)(MPI_Wtime() - iterstart);
#else
            double duration = (std::chrono::system_clock::now() - iterstart).count() / CLOCKS_PER_SEC;
#endif
            printiter(iter, drho, duration, diag_ethr);
            if (this->conv_elec)
            {
                this->niter = iter;
                bool stop = this->do_after_converge(iter);
                if (stop)
                    break;
            }
        }
        afterscf(istep);

        ModuleBase::timer::tick(this->classname, "Run");
    }

    return;
};

void ESolver_KS_LCAO_TDDFT::Run(const int istep, UnitCell_pseudo& ucell)
{
    if (!(GlobalV::CALCULATION == "scf" || GlobalV::CALCULATION == "md" || GlobalV::CALCULATION == "relax"
          || GlobalV::CALCULATION == "cell-relax" || GlobalV::CALCULATION.substr(0, 3) == "sto"))
    {
        this->othercalculation(istep);
    }
    else
    {
        ModuleBase::timer::tick(this->classname, "Run");

        if (this->maxniter > 0)
            this->printhead(); // print the headline on the screen.
        this->beforescf(istep); // Something else to do before the iter loop

        bool firstscf = true;
        this->conv_elec = false;
        this->niter = this->maxniter;
        for (int iter = 1; iter <= this->maxniter; ++iter)
        {
            writehead(GlobalV::ofs_running, istep, iter);
#ifdef __MPI
            auto iterstart = MPI_Wtime();
#else
            auto iterstart = std::chrono::system_clock::now();
#endif
            double diag_ethr = this->phsol->set_diagethr(istep, iter, drho);
            eachiterinit(istep, iter);
            this->hamilt2density(istep, iter, diag_ethr);

            //<Temporary> It may be changed when more clever parallel algorithm is put forward.
            // When parallel algorithm for bands are adopted. Density will only be treated in the first group.
            //(Different ranks should have abtained the same, but small differences always exist in practice.)
            // Maybe in the future, density and wavefunctions should use different parallel algorithms, in which
            // they do not occupy all processors, for example wavefunctions uses 20 processors while density uses 10.
            if (GlobalV::MY_STOGROUP == 0)
            {
                // double drho = this->estate.caldr2();
                // EState should be used after it is constructed.
                drho = GlobalC::CHR.get_drho();
                double hsolver_error = 0.0;
                if (firstscf)
                {
                    firstscf = false;
                    hsolver_error = this->phsol->cal_hsolerror();
                    // The error of HSolver is larger than drho, so a more precise HSolver should be excuconv_elected.
                    if (hsolver_error > drho)
                    {
                        diag_ethr = this->phsol->reset_diagethr(GlobalV::ofs_running, hsolver_error, drho);
                        this->hamilt2density(istep, iter, diag_ethr);
                        drho = GlobalC::CHR.get_drho();
                        hsolver_error = this->phsol->cal_hsolerror();
                    }
                }

                this->conv_elec = (drho < this->scf_thr);

                // If drho < hsolver_error in the first iter or drho < scf_thr, we do not change rho.
                if (drho < hsolver_error || this->conv_elec)
                {
                    if (drho < hsolver_error)
                        GlobalV::ofs_warning << " drho < hsolver_error, keep charge density unchanged." << std::endl;
                }
                else
                {
                    // charge mixing
                    // conv_elec = this->estate.mix_rho();
                    GlobalC::CHR.mix_rho(iter);
                }
            }
#ifdef __MPI
            MPI_Bcast(&drho, 1, MPI_DOUBLE, 0, PARAPW_WORLD);
            MPI_Bcast(&this->conv_elec, 1, MPI_DOUBLE, 0, PARAPW_WORLD);
            MPI_Bcast(GlobalC::CHR.rho[0], GlobalC::rhopw->nrxx, MPI_DOUBLE, 0, PARAPW_WORLD);
#endif

            // Hamilt should be used after it is constructed.
            // this->phamilt->update(conv_elec);
            updatepot(istep, iter);
            eachiterfinish(iter);
#ifdef __MPI
            double duration = (double)(MPI_Wtime() - iterstart);
#else
            double duration = (std::chrono::system_clock::now() - iterstart).count() / CLOCKS_PER_SEC;
#endif
            printiter(iter, drho, duration, diag_ethr);
            if (this->conv_elec)
            {
                this->niter = iter;
                bool stop = this->do_after_converge(iter);
                if (stop)
                    break;
            }
        }
        afterscf(istep);

        ModuleBase::timer::tick(this->classname, "Run");
    }

    return;
};

void ESolver_KS_LCAO_TDDFT::cal_Force(ModuleBase::matrix& force, ModuleBase::Vector3<double>* vel)
{
    Force_Stress_LCAO_TDDFT FSL(this->RA);
    FSL.getForceStress(GlobalV::CAL_FORCE,
                       GlobalV::CAL_STRESS,
                       GlobalV::TEST_FORCE,
                       GlobalV::TEST_STRESS,
                       this->LOC,
                       this->psid,
                       this->psi,
                       this->UHM,
                       force,
                       this->scs,
                       vel);
    // delete RA after cal_Force
    this->RA.delete_grid();
    this->have_force = true;
}

void ESolver_KS_LCAO_TDDFT::cal_Force(ModuleBase::matrix& force)
{
    Force_Stress_LCAO FSL(this->RA);
    FSL.getForceStress(GlobalV::CAL_FORCE,
                       GlobalV::CAL_STRESS,
                       GlobalV::TEST_FORCE,
                       GlobalV::TEST_STRESS,
                       this->LOC,
                       this->psid,
                       this->psi,
                       this->UHM,
                       force,
                       this->scs);
    // delete RA after cal_Force
    this->RA.delete_grid();
    this->have_force = true;
}

// use the original formula (Hamiltonian matrix) to calculate energy density matrix
void ESolver_KS_LCAO_TDDFT::cal_edm_tddft()
{
    this->LOC.edm_k_tddft.resize(GlobalC::kv.nks);
    for (int ik = 0; ik < GlobalC::kv.nks; ++ik)
    {
#ifdef __MPI
        this->LOC.edm_k_tddft[ik].create(this->LOC.ParaV->ncol, this->LOC.ParaV->nrow);
        complex<double>* Htmp = new complex<double>[this->LOC.ParaV->nloc];
        complex<double>* Sinv = new complex<double>[this->LOC.ParaV->nloc];
        complex<double>* tmp1 = new complex<double>[this->LOC.ParaV->nloc];
        complex<double>* tmp2 = new complex<double>[this->LOC.ParaV->nloc];
        complex<double>* tmp3 = new complex<double>[this->LOC.ParaV->nloc];
        complex<double>* tmp4 = new complex<double>[this->LOC.ParaV->nloc];
        complex<double>* tmp5 = new complex<double>[this->LOC.ParaV->nloc];
        ModuleBase::GlobalFunc::ZEROS(Htmp, this->LOC.ParaV->nloc);
        ModuleBase::GlobalFunc::ZEROS(Sinv, this->LOC.ParaV->nloc);
        ModuleBase::GlobalFunc::ZEROS(tmp1, this->LOC.ParaV->nloc);
        ModuleBase::GlobalFunc::ZEROS(tmp2, this->LOC.ParaV->nloc);
        ModuleBase::GlobalFunc::ZEROS(tmp3, this->LOC.ParaV->nloc);
        ModuleBase::GlobalFunc::ZEROS(tmp4, this->LOC.ParaV->nloc);
        ModuleBase::GlobalFunc::ZEROS(tmp5, this->LOC.ParaV->nloc);
        const int inc = 1;
        int nrow = this->LOC.ParaV->nrow;
        int ncol = this->LOC.ParaV->ncol;
        hamilt::MatrixBlock<complex<double>> h_mat, s_mat;
        p_hamilt->matrix(h_mat, s_mat);
        zcopy_(&this->LOC.ParaV->nloc, h_mat.p, &inc, Htmp, &inc);
        zcopy_(&this->LOC.ParaV->nloc, s_mat.p, &inc, Sinv, &inc);

        int* ipiv = new int[this->LOC.ParaV->nloc];
        int info;
        const int one_int = 1;
        pzgetrf_(&GlobalV::NLOCAL, &GlobalV::NLOCAL, Sinv, &one_int, &one_int, this->LOC.ParaV->desc, ipiv, &info);

        int LWORK = -1, liWORK = -1;
        std::vector<std::complex<double>> WORK(1, 0);
        std::vector<int> iWORK(1, 0);

        pzgetri_(&GlobalV::NLOCAL,
                 Sinv,
                 &one_int,
                 &one_int,
                 this->LOC.ParaV->desc,
                 ipiv,
                 WORK.data(),
                 &LWORK,
                 iWORK.data(),
                 &liWORK,
                 &info);

        LWORK = WORK[0].real();
        WORK.resize(LWORK, 0);
        liWORK = iWORK[0];
        iWORK.resize(liWORK, 0);

        pzgetri_(&GlobalV::NLOCAL,
                 Sinv,
                 &one_int,
                 &one_int,
                 this->LOC.ParaV->desc,
                 ipiv,
                 WORK.data(),
                 &LWORK,
                 iWORK.data(),
                 &liWORK,
                 &info);

        const char N_char = 'N', T_char = 'T';
        const double one_float[2] = {1.0, 0.0}, zero_float[2] = {0.0, 0.0};
        const complex<double> half_float[2] = {0.5, 0.0};
        pzgemm_(&T_char,
                &T_char,
                &GlobalV::NLOCAL,
                &GlobalV::NLOCAL,
                &GlobalV::NLOCAL,
                &one_float[0],
                this->LOC.dm_k[ik].c,
                &one_int,
                &one_int,
                this->LOC.ParaV->desc,
                Htmp,
                &one_int,
                &one_int,
                this->LOC.ParaV->desc,
                &zero_float[0],
                tmp1,
                &one_int,
                &one_int,
                this->LOC.ParaV->desc);

        pzgemm_(&N_char,
                &N_char,
                &GlobalV::NLOCAL,
                &GlobalV::NLOCAL,
                &GlobalV::NLOCAL,
                &one_float[0],
                tmp1,
                &one_int,
                &one_int,
                this->LOC.ParaV->desc,
                Sinv,
                &one_int,
                &one_int,
                this->LOC.ParaV->desc,
                &zero_float[0],
                tmp2,
                &one_int,
                &one_int,
                this->LOC.ParaV->desc);

        pzgemm_(&N_char,
                &T_char,
                &GlobalV::NLOCAL,
                &GlobalV::NLOCAL,
                &GlobalV::NLOCAL,
                &one_float[0],
                Sinv,
                &one_int,
                &one_int,
                this->LOC.ParaV->desc,
                Htmp,
                &one_int,
                &one_int,
                this->LOC.ParaV->desc,
                &zero_float[0],
                tmp3,
                &one_int,
                &one_int,
                this->LOC.ParaV->desc);

        pzgemm_(&N_char,
                &T_char,
                &GlobalV::NLOCAL,
                &GlobalV::NLOCAL,
                &GlobalV::NLOCAL,
                &one_float[0],
                tmp3,
                &one_int,
                &one_int,
                this->LOC.ParaV->desc,
                this->LOC.dm_k[ik].c,
                &one_int,
                &one_int,
                this->LOC.ParaV->desc,
                &zero_float[0],
                tmp4,
                &one_int,
                &one_int,
                this->LOC.ParaV->desc);

        pzgeadd_(&N_char,
                 &GlobalV::NLOCAL,
                 &GlobalV::NLOCAL,
                 &half_float[0],
                 tmp2,
                 &one_int,
                 &one_int,
                 this->LOC.ParaV->desc,
                 &half_float[0],
                 tmp4,
                 &one_int,
                 &one_int,
                 this->LOC.ParaV->desc);

        pztranu_(&GlobalV::NLOCAL,
                 &GlobalV::NLOCAL,
                 &one_float[0],
                 tmp4,
                 &one_int,
                 &one_int,
                 this->LOC.ParaV->desc,
                 &zero_float[0],
                 tmp5,
                 &one_int,
                 &one_int,
                 this->LOC.ParaV->desc);
        zcopy_(&this->LOC.ParaV->nloc, tmp5, &inc, this->LOC.edm_k_tddft[ik].c, &inc);

        delete[] Htmp;
        delete[] Sinv;
        delete[] tmp1;
        delete[] tmp2;
        delete[] tmp3;
        delete[] tmp4;
        delete[] tmp5;
        delete[] ipiv;
#else
        this->LOC.edm_k_tddft[ik].create(this->LOC.ParaV->ncol, this->LOC.ParaV->nrow);
        ModuleBase::ComplexMatrix Sinv(GlobalV::NLOCAL, GlobalV::NLOCAL);
        ModuleBase::ComplexMatrix Htmp(GlobalV::NLOCAL, GlobalV::NLOCAL);
        hamilt::MatrixBlock<complex<double>> h_mat, s_mat;
        p_hamilt->matrix(h_mat, s_mat);
        // cout<<"hmat "<<h_mat.p[0]<<endl;
        for (int i = 0; i < GlobalV::NLOCAL; i++)
        {
            for (int j = 0; j < GlobalV::NLOCAL; j++)
            {
                Htmp(i, j) = h_mat.p[i * GlobalV::NLOCAL + j];
                Sinv(i, j) = s_mat.p[i * GlobalV::NLOCAL + j];
            }
        }
        int INFO;

        int LWORK = 3 * GlobalV::NLOCAL - 1; // tmp
        std::complex<double>* WORK = new std::complex<double>[LWORK];
        ModuleBase::GlobalFunc::ZEROS(WORK, LWORK);
        int IPIV[GlobalV::NLOCAL];

        LapackConnector::zgetrf(GlobalV::NLOCAL, GlobalV::NLOCAL, Sinv, GlobalV::NLOCAL, IPIV, &INFO);
        LapackConnector::zgetri(GlobalV::NLOCAL, Sinv, GlobalV::NLOCAL, IPIV, WORK, LWORK, &INFO);

        this->LOC.edm_k_tddft[ik] = 0.5 * (Sinv * Htmp * this->LOC.dm_k[ik] + this->LOC.dm_k[ik] * Htmp * Sinv);
        delete[] WORK;
#endif
    }
    return;
}
} // namespace ModuleESolver