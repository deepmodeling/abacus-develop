#include "esolver_of.h"

#include "module_io/rho_io.h"
#include "module_io/potential_io.h"
#include "module_io/output_log.h"
//-----------temporary-------------------------
#include "module_base/global_function.h"
#include "module_base/memory.h"
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_hamilt_general/module_ewald/H_Ewald_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/print_info.h"
//-----force-------------------
#include "module_hamilt_pw/hamilt_pwdft/forces.h"
//-----stress------------------
#include "module_hamilt_pw/hamilt_ofdft/of_stress_pw.h"
//---------------------------------------------------
#include "module_hamilt_pw/hamilt_pwdft/hamilt_pw.h"

namespace ModuleESolver
{

void ESolver_OF::Init(Input &inp, UnitCell &ucell)
{
    ESolver_FP::Init(inp, ucell);

    // save necessary parameters
    this->of_kinetic_ = inp.of_kinetic;
    this->of_method_ = inp.of_method;
    this->of_conv_ = inp.of_conv;
    this->of_tole_ = inp.of_tole;
    this->of_tolp_ = inp.of_tolp;
    this->max_iter_ = inp.scf_nmax;

    ucell.cal_nelec(GlobalV::nelec);

	if(ucell.atoms[0].ncpp.xc_func=="HSE"||ucell.atoms[0].ncpp.xc_func=="PBE0")
	{
        ModuleBase::WARNING_QUIT("esolver_of", "Hybrid functionals are not supported by OFDFT.");
		// XC_Functional::set_xc_type("pbe");
	}
	else
	{
		XC_Functional::set_xc_type(ucell.atoms[0].ncpp.xc_func);
	}

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SETUP UNITCELL");

    // symmetry analysis should be performed every time the cell is changed
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        this->symm.analy_sys(ucell, GlobalV::ofs_running);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SYMMETRY");
    }

    // Setup the k points according to symmetry.
    kv.set(this->symm, GlobalV::global_kpoint_card, GlobalV::NSPIN, ucell.G, ucell.latvec);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"INIT K-POINTS");

    // print information
    // mohan add 2021-01-30
    Print_Info::setup_parameters(ucell, kv);

    // initialize the real-space uniform grid for FFT and parallel
    // distribution of plane waves
    GlobalC::Pgrid.init(pw_rho->nx,
                        pw_rho->ny,
                        pw_rho->nz,
                        pw_rho->nplane,
                        pw_rho->nrxx,
                        pw_big->nbz,
                        pw_big->bz); // mohan add 2010-07-22, update 2011-05-04
    // Calculate Structure factor
    sf.setup_structure_factor(&GlobalC::ucell, pw_rho);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT BASIS");

    this->dV_ = ucell.omega / this->pw_rho->nxyz; // volume of one point in real space

    //----------------------------------------------------------
    // 1 read in initial data:
    //   a lattice structure:atom_species,atom_positions,lattice vector
    //   b k_points
    //   c pseudopotential
    // 2 setup planeware basis, FFT,structure factor, ...
    // 3 initialize local pseudopotential in G_space
    // 4 initialize charge desity and warefunctios in real space
    //----------------------------------------------------------

    // Initialize the "wavefunction", which is sqrt(rho)
    this->psi_ = new psi::Psi<double>(1, GlobalV::NSPIN, this->pw_rho->nrxx);
    ModuleBase::Memory::record("OFDFT::Psi", sizeof(double) * GlobalV::NSPIN * this->pw_rho->nrxx);
    this->pphi_ = new double*[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        this->pphi_[is] = this->psi_->get_pointer(is);
    }
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT PHI");


    //=================================
    // initalize local pseudopotential
    //=================================
    GlobalC::ppcell.init_vloc(GlobalC::ppcell.vloc, pw_rho);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "LOCAL POTENTIAL");

    //======================================
    // Initalize non local pseudopotential
    //======================================
    GlobalC::ppcell.init_vnl(GlobalC::ucell, pw_rho);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "NON-LOCAL POTENTIAL");

    if(this->pelec == nullptr)
    {
        this->pelec = new elecstate::ElecState((Charge*)(&chr), this->pw_rho, pw_big);
    }

    this->pelec->charge->allocate(GlobalV::NSPIN);
    this->pelec->omega = GlobalC::ucell.omega;

    this->pelec->pot = new elecstate::Potential(pw_rhod,
                                                pw_rho,
                                                &GlobalC::ucell,
                                                &GlobalC::ppcell.vloc,
                                                &sf,
                                                &(this->pelec->f_en.etxc),
                                                &(this->pelec->f_en.vtxc));
    //There is no Operator in ESolver_OF, register Potentials here!
    std::vector<std::string> pot_register_in;
    if (GlobalV::VION_IN_H)
    {
        pot_register_in.push_back("local");
    }
    if (GlobalV::VH_IN_H)
    {
        pot_register_in.push_back("hartree");
    }
    //no variable can choose xc, maybe it is necessary
    pot_register_in.push_back("xc");
    if (GlobalV::imp_sol)
    {
        pot_register_in.push_back("surchem");
    }
    if (GlobalV::EFIELD_FLAG)
    {
        pot_register_in.push_back("efield");
    }
    if (GlobalV::GATE_FLAG)
    {
        pot_register_in.push_back("gatefield");
    }
    //only Potential is not empty, Veff and Meta are available
    if(pot_register_in.size()>0)
    {
        //register Potential by gathered operator
        this->pelec->pot->pot_register(pot_register_in);
    }

    //=========================================================
    // calculate the total local pseudopotential in real space
    //=========================================================
    this->pelec->init_scf(0, sf.strucFac); // atomic_rho, v_of_rho, set_vrs

    // liuyu move here 2023-10-09
    // D in uspp need vloc, thus behind init_scf()
    // calculate the effective coefficient matrix for non-local pseudopotential projectors
    ModuleBase::matrix veff = this->pelec->pot->get_effective_v();
    GlobalC::ppcell.cal_effective_D(veff, this->pw_rho, GlobalC::ucell);

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT POTENTIAL");

    // Calculate electron numbers
    this->nelec_ = new double[GlobalV::NSPIN];
    if (GlobalV::NSPIN == 1)
    {
        this->nelec_[0] = GlobalV::nelec;
    }
    else if (GlobalV::NSPIN == 2)
    {
        //in fact, nelec_spin will not be used anymore
        this->pelec->init_nelec_spin();
        this->nelec_[0] = this->pelec->nelec_spin[0];
        this->nelec_[1] = this->pelec->nelec_spin[1];
    }

    // ================================
    // Initialize optimization methods
    // ================================
    this->init_opt();

    // optimize theta if nspin=2
    if (GlobalV::NSPIN == 2)
    {
        this->opt_cg_mag_ = new ModuleBase::Opt_CG;
        this->opt_cg_mag_->allocate(GlobalV::NSPIN);
    }

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT OPTIMIZATION");

    // =============================================
    // Initalize chemical potential, step length, ...
    // =============================================
    this->mu_ = new double[GlobalV::NSPIN];
    this->theta_ = new double[GlobalV::NSPIN];
    this->pdLdphi_ = new double*[GlobalV::NSPIN];
    this->pdEdphi_ = new double*[GlobalV::NSPIN];
    this->pdirect_ = new double*[GlobalV::NSPIN];
    this->precipDir_ = new std::complex<double> *[GlobalV::NSPIN];
    // this->pdeltaRhoHar = new double[this->pw_rho->nrxx];

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        this->pdLdphi_[is] = new double[this->pw_rho->nrxx];
        this->pdEdphi_[is] = new double[this->pw_rho->nrxx];
        this->pdirect_[is] = new double[this->pw_rho->nrxx];
        this->precipDir_[is] = new std::complex<double>[pw_rho->npw];
    }

    // ===================================
    // Initialize KEDF
    // ===================================
    this->init_kedf();

    // Initialize charge extrapolation
    CE.Init_CE(GlobalC::ucell.nat);
    delete this->ptempRho_;
    this->ptempRho_ = new Charge();
    this->ptempRho_->set_rhopw(this->pw_rho);
    this->ptempRho_->allocate(GlobalV::NSPIN);
}

void ESolver_OF::init_after_vc(Input &inp, UnitCell &ucell)
{
    ModuleBase::timer::tick("ESolver_OF", "init_after_vc");

    ESolver_FP::init_after_vc(inp,ucell);

    this->pw_rho->nrxx = this->pw_rho->nrxx;
    this->dV_ = ucell.omega / this->pw_rho->nxyz; // volume of one point in real space

    if (GlobalV::md_prec_level == 2)
    {
        // initialize the real-space uniform grid for FFT and parallel
        // distribution of plane waves
        GlobalC::Pgrid.init(this->pw_rho->nx,
                            this->pw_rho->ny,
                            this->pw_rho->nz,
                            this->pw_rho->nplane,
                            this->pw_rho->nrxx,
                            pw_big->nbz,
                            pw_big->bz); // mohan add 2010-07-22, update 2011-05-04

        // Calculate Structure factor
        this->sf.setup_structure_factor(&ucell, this->pw_rho);
    }

    delete this->pelec;
    this->pelec = new elecstate::ElecState((Charge*)(&chr), this->pw_rho, pw_big);

    this->pelec->charge->allocate(GlobalV::NSPIN);
    this->pelec->omega = GlobalC::ucell.omega;

    delete this->pelec->pot;
    this->pelec->pot = new elecstate::Potential(this->pw_rhod,
                                                this->pw_rho,
                                                &GlobalC::ucell,
                                                &(GlobalC::ppcell.vloc),
                                                &(this->sf),
                                                &(this->pelec->f_en.etxc),
                                                &(this->pelec->f_en.vtxc));
    //There is no Operator in ESolver_OF, register Potentials here!
    std::vector<std::string> pot_register_in;
    if (GlobalV::VION_IN_H)
    {
        pot_register_in.push_back("local");
    }
    if (GlobalV::VH_IN_H)
    {
        pot_register_in.push_back("hartree");
    }
    //no variable can choose xc, maybe it is necessary
    pot_register_in.push_back("xc");
    if (GlobalV::imp_sol)
    {
        pot_register_in.push_back("surchem");
    }
    if (GlobalV::EFIELD_FLAG)
    {
        pot_register_in.push_back("efield");
    }
    if (GlobalV::GATE_FLAG)
    {
        pot_register_in.push_back("gatefield");
    }
    //only Potential is not empty, Veff and Meta are available
    if(pot_register_in.size()>0)
    {
        //register Potential by gathered operator
        this->pelec->pot->pot_register(pot_register_in);
    }

    // ================================
    // Initialize optimization methods
    // ================================
    this->init_opt();

    GlobalC::ppcell.init_vnl(GlobalC::ucell, pw_rho);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"NON-LOCAL POTENTIAL");

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        delete[] this->pdLdphi_[is];
        delete[] this->pdEdphi_[is];
        delete[] this->pdirect_[is];
        delete[] this->precipDir_[is];
        this->pdLdphi_[is] = new double[this->pw_rho->nrxx];
        this->pdEdphi_[is] = new double[this->pw_rho->nrxx];
        this->pdirect_[is] = new double[this->pw_rho->nrxx];
        this->precipDir_[is] = new std::complex<double>[pw_rho->npw];
    }

    // ===================================
    // Initialize KEDF
    // ===================================
    this->init_kedf();

    delete this->ptempRho_;
    this->ptempRho_ = new Charge();
    this->ptempRho_->set_rhopw(this->pw_rho);
    this->ptempRho_->allocate(GlobalV::NSPIN);
}

void ESolver_OF::Run(int istep, UnitCell& ucell)
{
    ModuleBase::timer::tick("ESolver_OF", "Run");
    // get Ewald energy, initial rho and phi if necessary
    this->beforeOpt(istep);
    this->iter_ = 0;

    while(true)
    {
        // once we get a new rho and phi, update potential
        this->updateV();

        // calculate the energy of new rho and phi
        this->energy_llast_ = this->energy_last_;
        this->energy_last_ = this->energy_current_;
        this->energy_current_ = this->cal_Energy();

        // print neccesary information
        this->printInfo();

        // check if the job is done
        if (this->checkExit()) break;

        // find the optimization direction and step lenghth theta according to the potential
        this->solveV();

        // update the rho and phi based on the direction and theta
        this->updateRho();

        this->iter_++;
    }

    this->afterOpt(istep);

    ModuleBase::timer::tick("ESolver_OF", "Run");
}

//
// Calculate ewald energy, initialize the rho, phi, theta
//
void ESolver_OF::beforeOpt(const int istep)
{
    if (GlobalC::ucell.cell_parameter_updated)
    {
        this->init_after_vc(INPUT, GlobalC::ucell);
    }
    if (GlobalC::ucell.ionic_position_updated)
    {
        CE.update_all_dis(GlobalC::ucell);
        CE.extrapolate_charge(
#ifdef __MPI
            &(GlobalC::Pgrid),
#endif
            GlobalC::ucell,
            pelec->charge,
            &(sf));
    }

    this->pelec->init_scf(istep, sf.strucFac);

    //calculate ewald energy
    this->pelec->f_en.ewald_energy = H_Ewald_pw::compute_ewald(GlobalC::ucell, this->pw_rho, sf.strucFac);

    Symmetry_rho srho;
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        srho.begin(is, *(pelec->charge), this->pw_rho, GlobalC::Pgrid, this->symm);
    }

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        if (GlobalV::init_chg != "file")
        {
            for (int ibs = 0; ibs < this->pw_rho->nrxx; ++ibs)
            {
                // Here we initialize rho to be uniform,
                // because the rho got by pot.init_pot -> Charge::atomic_rho may contain minus elements.
                pelec->charge->rho[is][ibs] = this->nelec_[is]/GlobalC::ucell.omega;
                this->pphi_[is][ibs] = sqrt(pelec->charge->rho[is][ibs]);
            }
        }
        else
        {
            for (int ibs = 0; ibs < this->pw_rho->nrxx; ++ibs)
            {
                this->pphi_[is][ibs] = sqrt(pelec->charge->rho[is][ibs]);
            }
        }
    }

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        this->mu_[is] = 0;
        this->theta_[is] = 0.;
        ModuleBase::GlobalFunc::ZEROS(this->pdLdphi_[is], this->pw_rho->nrxx);
        ModuleBase::GlobalFunc::ZEROS(this->pdEdphi_[is], this->pw_rho->nrxx);
        ModuleBase::GlobalFunc::ZEROS(this->pdirect_[is], this->pw_rho->nrxx);
    }
    if (GlobalV::NSPIN == 1)
    {
        this->theta_[0] = 0.2;
    }
}

//
// Get dL/dphi = dL/drho * drho/dphi = (dE/drho - mu) * 2 * phi
//
void ESolver_OF::updateV()
{
    // (1) get dL/dphi
    if(GlobalV::NSPIN==4) GlobalC::ucell.cal_ux();
    this->pelec->pot->update_from_charge(pelec->charge, &GlobalC::ucell); // Hartree + XC + external
    this->kineticPotential(pelec->charge->rho, this->pphi_, this->pelec->pot->get_effective_v()); // (kinetic + Hartree + XC + external) * 2 * phi
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        const double* vr_eff = this->pelec->pot->get_effective_v(is);
        for (int ir = 0; ir < this->pw_rho->nrxx; ++ir)
        {
            this->pdEdphi_[is][ir] = vr_eff[ir];
        }
        this->mu_[is] = this->cal_mu(this->pphi_[is], this->pdEdphi_[is], this->nelec_[is]);

        for (int ir = 0; ir < this->pw_rho->nrxx; ++ir)
        {
            this->pdLdphi_[is][ir] = this->pdEdphi_[is][ir] - 2. * this->mu_[is] * this->pphi_[is][ir];
        }
    }

    // (2) get the norm of dLdphi
    // ===== temporary solution of potential convergence when of_full_pw = 0 =====
    this->normdLdphi_llast_ = this->normdLdphi_last_;
    this->normdLdphi_last_ = this->normdLdphi_;
    // ===========================================================================
    this->normdLdphi_ = 0.;

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
       this->normdLdphi_ += this->inner_product(this->pdLdphi_[is], this->pdLdphi_[is], this->pw_rho->nrxx, 1);
    }
    Parallel_Reduce::reduce_all(this->normdLdphi_);
    this->normdLdphi_ = sqrt(this->normdLdphi_/this->pw_rho->nxyz/GlobalV::NSPIN);
}

//
// Get optimization direction d and step theta
//
void ESolver_OF::solveV()
{
    // (1) get |d0> with optimization algorithm
    this->get_direction();
    // initialize tempPhi and tempRho used in line search
    double **ptempPhi = new double*[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        ptempPhi[is] = new double[this->pw_rho->nrxx];
        for (int ir = 0; ir < this->pw_rho->nrxx; ++ir)
        {
            ptempPhi[is][ir] = this->pphi_[is][ir];
            this->ptempRho_->rho[is][ir] = ptempPhi[is][ir] * ptempPhi[is][ir];
        }
    }

    // (2) rotate and renormalize the direction
    this->getNextDirect();

    // (3) make sure dEdtheta<0 at theta = 0
    double E = 0.; // energy of tempPhi and tempRho
    double *dEdtheta = new double[GlobalV::NSPIN]; // dE/dtheta of tempPhi
    double *tempTheta = new double[GlobalV::NSPIN];
    ModuleBase::GlobalFunc::ZEROS(dEdtheta, GlobalV::NSPIN);
    ModuleBase::GlobalFunc::ZEROS(tempTheta, GlobalV::NSPIN);

    double dEdthetaThre = 1e5; // threshould of dEdtheta, avoid the unstable optimization
    this->caldEdtheta(ptempPhi, this->ptempRho_, tempTheta, dEdtheta);

    // Assert dEdtheta(theta = 0) < 0, otherwise line search will not work.
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        if (dEdtheta[is] > dEdthetaThre)
        {
            std::cout << "dEdtheta    " << dEdtheta[is] << std::endl;
            ModuleBase::WARNING_QUIT("esolver_of.cpp", "dE/dtheta is too large.");
        }
        else if (dEdtheta[is] > 0)
        {
            GlobalV::ofs_warning << "ESolver_OF: WARNING " << "dEdphi > 0, replace direct with steepest descent method." << std::endl;
            for (int ir = 0; ir < this->pw_rho->nrxx; ++ir)
            {
                this->pdirect_[is][ir] = - this->pdLdphi_[is][ir];
            }
            this->getNextDirect();
            this->caldEdtheta(ptempPhi, this->ptempRho_, tempTheta, dEdtheta);
            if (dEdtheta[is] > dEdthetaThre)
            {
                std::cout << "dEdtheta    " << dEdtheta[is] << std::endl;
                ModuleBase::WARNING_QUIT("esolver_of.cpp", "dE/dtheta is too large.");
            }
            else if (dEdtheta[is] > 0)
            {
                GlobalV::ofs_warning << "ESolver_OF: WARNING " << "when use steepest dencent method, dEdphi > 0, so we might get minimum." << std::endl;
            }
        }
    }
    delete[] tempTheta;

    // // ======================== for test ============================
    //     if (this->iter_ == 0)
    //     {
    //         for (int i = -100; i < 100; ++i)
    //         {
    //             this->theta_[0] = 0.001 * i;
    //             for (int ir = 0; ir < this->pw_rho->nrxx; ++ir)
    //             {
    //                 ptempPhi[0][ir] = this->pphi_[0][ir] * cos(this->theta_[0]) + this->pdirect_[0][ir] *
    //                 sin(this->theta_[0]); ptempRho->rho[0][ir] = ptempPhi[0][ir] * ptempPhi[0][ir];
    //             }
    //             this->caldEdtheta(ptempPhi, ptempRho, this->theta_, dEdtheta);
    //             this->pelec->f_en.calculate_etot(this->pw_rho->nrxx, this->pw_rho->nxyz);
    //             E = this->pelec->f_en.etot;
    //             double eKE = 0.;
    //             double ePP = 0.;
    //             eKE = this->kineticEnergy();
    //             ePP = this->inner_product(this->pelec->pot->get_fixed_v(), ptempRho->rho[0], this->pw_rho->nrxx, this->dV_);
    //             // ePP = this->inner_product(GlobalC::pot.vltot, ptempRho[0], this->pw_rho->nrxx, this->dV_);
    //             Parallel_Reduce::reduce_all(ePP);
    //             E += eKE + ePP;
    //             GlobalV::ofs_warning << i << "    " << dEdtheta[0] << "    " << E << endl;
    //             if (this->theta_[0] == 0) cout << "dEdtheta    " << dEdtheta[0]<< endl;
    //         }
    //         exit(0);
    //     }
    // // ======================== for test ============================

    // (4) line search to find the best theta
    double eKE = 0.;    // kinetic energy
    double ePP = 0.;    // electron-ion interaction energy
    if (GlobalV::NSPIN == 1)
    {
        int numDC = 0; // iteration number of line search
        strcpy(this->task_, "START");
        while (true)
        {
            // update energy
            this->pelec->cal_energies(2);
            E = this->pelec->f_en.etot;
            eKE = this->kineticEnergy();
            ePP = this->inner_product(this->pelec->pot->get_fixed_v(), this->ptempRho_->rho[0], this->pw_rho->nrxx, this->dV_);
            Parallel_Reduce::reduce_all(ePP);
            E += eKE + ePP;

            // line search to update theta[0]
            this->opt_dcsrch_->dcSrch(E, dEdtheta[0], this->theta_[0], this->task_);
            numDC++;

            // decide what to do next according to the output of line search
            if (strncmp(this->task_, "FG", 2) == 0) // continue line search
            {
                // update tempPhi and tempRho
                for (int i = 0; i < this->pw_rho->nrxx; ++i)
                {
                    ptempPhi[0][i] = this->pphi_[0][i] * cos(this->theta_[0]) + this->pdirect_[0][i] * sin(this->theta_[0]);
                    this->ptempRho_->rho[0][i] = ptempPhi[0][i] * ptempPhi[0][i];
                }
                // get dEdtheta of new tempPhi and tempRho
                this->caldEdtheta(ptempPhi, this->ptempRho_, this->theta_, dEdtheta);

                if (numDC > this->maxDCsrch_)
                {
                    GlobalV::ofs_warning << "ESolver_OF linesearch: WARNING " << "excedd the max iter number." << std::endl;
                    break;
                }
            }
            else if (strncmp(this->task_, "CO", 2) == 0) // convergence achieved
            {
                break;
            }
            else if (strncmp(this->task_, "WA", 2) == 0) // warning of line search
            {
                GlobalV::ofs_warning << "ESolver_OF linesearch: WARNING " << this->task_ << std::endl;
                std::cout << this->task_ << std::endl;
                break;
            }
            else if (strncmp(this->task_, "ER", 2) == 0) // ERROR in line search
            {
                GlobalV::ofs_warning << "ESolver_OF linesearch: ERROR " << this->task_ << std::endl;
                std::cout << this->task_ << std::endl;
                break;
            }
        }
    }
    else if (GlobalV::NSPIN == 2)
    {
        ModuleBase::WARNING_QUIT("esolver_of", "Sorry, SPIN2 case is not supported by OFDFT for now.");
    // ========================== Under testing ==========================
    //     this->opt_cg_mag_->refresh();

    //     double *pthetaDir = new double[GlobalV::NSPIN];
    //     double *tempTheta = new double[GlobalV::NSPIN];
    //     ModuleBase::GlobalFunc::ZEROS(pthetaDir, GlobalV::NSPIN);
    //     ModuleBase::GlobalFunc::ZEROS(tempTheta, GlobalV::NSPIN);
    //     double thetaAlpha = 0.;
    //     double alphaTol = 1e-4;
    //     double maxThetaDir = 0.;
    //     double dEdalpha = 0.;
    //     int thetaIter = 0;
    //     int numDC = 0;

    //     while (true)
    //     {
    //         this->opt_cg_mag_->next_direct(dEdtheta, 1, pthetaDir);

    //         dEdalpha = this->inner_product(dEdtheta, pthetaDir, 2, 1.);

    //         if (dEdalpha >= 0.)
    //         {
    //             for (int is = 0; is < GlobalV::NSPIN; ++is)
    //             {
    //                 pthetaDir[is] = -dEdtheta[is];
    //             }
    //             dEdalpha = this->inner_product(dEdtheta, pthetaDir, 2, 1);
    //         }

    //         maxThetaDir = max(abs(pthetaDir[0]), abs(pthetaDir[1]));
    //         thetaAlpha = min(0.1, 0.1*ModuleBase::PI/maxThetaDir);

        //         // line search along thetaDir to find thetaAlpha
        //         this->opt_dcsrch_->set_paras(1e-4, 1e-2, 1e-12, 0., ModuleBase::PI/maxThetaDir);
        //         strcpy(this->task_, "START");
        //         numDC = 0;
        //         while(true)
        //         {
        //             this->pelec->f_en.calculate_etot(this->pw_rho->nrxx, this->pw_rho->nxyz);
        //             E = this->pelec->f_en.etot;
        //             eKE = this->kineticEnergy();
        //             ePP = 0.;
        //             for (int is = 0; is < GlobalV::NSPIN; ++is) {
        //                 ePP += this->inner_product(GlobalC::pot.vltot, ptempRho[is], this->pw_rho->nrxx, this->dV_);
        //             }
        //             Parallel_Reduce::reduce_all(ePP);
        //             E += eKE + ePP;
        //             this->opt_dcsrch_->dcSrch(E, dEdalpha, thetaAlpha, this->task_);
        //             numDC++;

        //             if (strncmp(this->task_, "FG", 2) == 0)
        //             {
        //                 for (int is = 0; is < GlobalV::NSPIN; ++is)
        //                 {
        //                     tempTheta[is] = this->theta_[is] + thetaAlpha * pthetaDir[is];
        //                     for (int ir = 0; ir < this->pw_rho->nrxx; ++ir)
        //                     {
        //                         ptempPhi[is][ir] = this->pphi_[is][ir] * cos(tempTheta[is]) + this->pdirect_[is][ir] *
        //                         sin(tempTheta[is]); ptempRho[is][ir] = ptempPhi[is][ir] * ptempPhi[is][ir];
        //                     }
        //                 }
        //                 this->caldEdtheta(ptempPhi, ptempRho, tempTheta, dEdtheta);
        //                 dEdalpha = this->inner_product(dEdtheta, pthetaDir, 2, 1);

        //                 if (numDC > 10)
        //                 {
        //                     GlobalV::ofs_warning << "ESolver_OF linesearch: WARNING " << "excedd the max iter
        //                     number." << endl; break;
        //                 }
        //             }
        //             else if (strncmp(this->task_, "CO", 2) == 0)
        //             {
        //                 break;
        //             }
        //             else if (strncmp(this->task_, "WA", 2) == 0)
        //             {
        //                 GlobalV::ofs_warning << "ESolver_OF linesearch: WARNING " << this->task_ << std::endl;
        //                 cout << this->task_ << endl;
        //                 break;
        //             }
        //             else if (strncmp(this->task_, "ER", 2) == 0)
        //             {
        //                 GlobalV::ofs_warning << "ESolver_OF linesearch: ERROR " << this->task_ << std::endl;
        //                 cout << this->task_ << endl;
        //                 break;
        //             }
        //         }

        //         for (int is = 0; is < GlobalV::NSPIN; ++is) this->theta_[is] += thetaAlpha * pthetaDir[is];
        //         if (sqrt(dEdtheta[0] * dEdtheta[0] + dEdtheta[1] * dEdtheta[1]) < alphaTol) break;
        //         thetaIter++;
        //         if (thetaIter > 2) break;
        //     }
        //     delete[] tempTheta;
        //     delete[] pthetaDir;
        // ========================== Under testing ==========================
    }

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        delete[] ptempPhi[is];
    }
    delete[] ptempPhi;
    delete[] dEdtheta;
}

//
// Rotate and renormalize the direction |d>, make it orthogonal to phi, and <d|d> = nelec
//
void ESolver_OF::getNextDirect()
{
    // filter the high frequency term in direction if of_full_pw = false
    if (!GlobalV::of_full_pw)
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            pw_rho->real2recip(this->pdirect_[is], this->precipDir_[is]);
            pw_rho->recip2real(this->precipDir_[is], this->pdirect_[is]);
        }
    }

    if (GlobalV::NSPIN == 1)
    {
        double tempTheta = 0; // tempTheta = |d'|/|d0 + phi|, theta = min(theta, tempTheta)

        // (1) make direction orthogonal to phi
        // |d'> = |d0> - |phi><phi|d0>/nelec
        double innerPhiDir = this->inner_product(this->pdirect_[0], this->pphi_[0], this->pw_rho->nrxx, this->dV_);
        Parallel_Reduce::reduce_all(innerPhiDir);
        for (int i = 0; i < this->pw_rho->nrxx; ++i)
        {
            tempTheta += pow(this->pdirect_[0][i] + this->pphi_[0][i], 2);
            this->pdirect_[0][i] = this->pdirect_[0][i] - this->pphi_[0][i] * innerPhiDir / this->nelec_[0];
        }
        Parallel_Reduce::reduce_all(tempTheta);
        tempTheta = sqrt(tempTheta);

        // (2) renormalize direction
        // |d> = |d'> * \sqrt(nelec) / <d'|d'>
        double normDir = this->inner_product(this->pdirect_[0], this->pdirect_[0], this->pw_rho->nrxx, this->dV_);
        Parallel_Reduce::reduce_all(normDir);
        normDir = sqrt(normDir);
        for (int i = 0; i < this->pw_rho->nrxx; ++i)
        {
            this->pdirect_[0][i] = sqrt(this->nelec_[0]) * this->pdirect_[0][i] / normDir;
        }

        tempTheta = normDir/tempTheta;
        this->theta_[0] = std::min(this->theta_[0], tempTheta);
    }
    else if (GlobalV::NSPIN == 2) // theta = 0
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            // (1) make direction orthogonal to phi
            // |d'> = |d0> - |phi><phi|d0>/nelec
            double innerPhiDir = this->inner_product(this->pdirect_[is], this->pphi_[is], this->pw_rho->nrxx, this->dV_);
            Parallel_Reduce::reduce_all(innerPhiDir);
            for (int i = 0; i < this->pw_rho->nrxx; ++i)
            {
                this->pdirect_[is][i] = this->pdirect_[is][i] - this->pphi_[is][i] * innerPhiDir / this->nelec_[is];
            }

            // (2) renormalize direction
            // |d> = |d'> * \sqrt(nelec) / <d'|d'>
            double normDir = this->inner_product(this->pdirect_[is], this->pdirect_[is], this->pw_rho->nrxx, this->dV_);
            Parallel_Reduce::reduce_all(normDir);
            normDir = sqrt(normDir);
            for (int i = 0; i < this->pw_rho->nrxx; ++i)
            {
                this->pdirect_[is][i] = sqrt(this->nelec_[is]) * this->pdirect_[is][i] / normDir;
            }
            this->theta_[is] = 0.;
        }
    }
}

//
// Update the density and "wavefunction" after one step of optimization
//
void ESolver_OF::updateRho()
{
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        for (int ir = 0; ir < this->pw_rho->nrxx; ++ir)
        {
            this->pphi_[is][ir] = this->pphi_[is][ir] * cos(this->theta_[is]) + this->pdirect_[is][ir] * sin(this->theta_[is]);
            pelec->charge->rho[is][ir] = this->pphi_[is][ir] * this->pphi_[is][ir];
        }
    }
    // ============================ for test ===========================
    // if (ModuleSymmetry::Symmetry::symm_flag == 1)
    // {
    //     Symmetry_rho srho;
    //     for (int is = 0; is < GlobalV::NSPIN; is++)
    //     {
    //         srho.begin(is, pelec->charge, this->pw_rho, GlobalC::Pgrid, this->symm);
    //         for (int ibs = 0; ibs < this->pw_rho->nrxx; ++ibs)
    //         {
    //             this->pphi_[is][ibs] = sqrt(pelec->charge->rho[is][ibs]);
    //         }
    //     }
    // }
    // =============test for rho convergence criterion =================
    // for (int is = 0; is < GlobalV::NSPIN; ++is)
    // {
    //     for (int ir = 0; ir < this->pw_rho->nrxx; ++ir)
    //     {
    //         pelec->charge->rho_save[is][ir] = pelec->charge->rho[is][ir];
    //         this->pdeltaRho[is][ir] = pelec->charge->rho[is][ir] - pelec->charge->rho_save[is][ir];
    //         this->deltaRhoR += abs(this->pdeltaRho[is][ir]);

    //     }
    //     this->pw_rho->real2recip(this->pdeltaRho[is], this->precipDir_[is]);
    //     for (int ig = 0; ig < this->pw_rho->npw; ++ig)
    //     {
    //         if (this->pw_rho->gg[ig] != 0.)
    //             this->precipDir_[is][ig] = this->precipDir_[is][ig] / this->pw_rho->gg[ig] / this->pw_rho->tpiba2 * 4. * M_PI;
    //         else
    //             this->precipDir_[is][ig] = 0.;
    //     }
    //     this->pw_rho->recip2real(this->precipDir_[is], this->pdeltaRhoHar);
    //     this->deltaRhoG = this->inner_product(this->pdeltaRho[is], this->pdeltaRhoHar, this->pw_rho->nrxx, this->dV_);
    // }
    // Parallel_Reduce::reduce_all(this->deltaRhoR);
    // Parallel_Reduce::reduce_all(this->deltaRhoG);
    // this->deltaRhoR *= this->dV_;
    // this->deltaRhoG /= 2.;
}

//
// Check convergence, return ture if converge or iter >= maxIter.
//
bool ESolver_OF::checkExit()
{
    this->conv_ = false;
    bool potConv = false;
    bool potHold = false; // if normdLdphi nearly remains unchanged
    bool energyConv = false;

    if (this->normdLdphi_ < this->of_tolp_)
        potConv = true;
    if (this->iter_ >= 3
        && std::abs(this->normdLdphi_ - this->normdLdphi_last_) < 1e-10
        && std::abs(this->normdLdphi_ - this->normdLdphi_llast_) < 1e-10)
        potHold = true;

    if (this->iter_ >= 3
        && std::abs(this->energy_current_ - this->energy_last_) < this->of_tole_
        && std::abs(this->energy_current_ - this->energy_llast_) < this->of_tole_)
        energyConv = true;

    if (this->of_conv_ == "energy" && energyConv)
    {
        this->conv_ = true;
        return true;
    }
    else if (this->of_conv_ == "potential" && potConv)
    {
        this->conv_ = true;
        return true;
    // ============ temporary solution of potential convergence ===========
    }
    else if (this->of_conv_ == "potential" && potHold)
    {
        GlobalV::ofs_warning << "ESolver_OF WARNING: " <<
        "The convergence of potential has not been reached, but the norm of potential nearly remains unchanged, set of_full_pw = 1 may work." << std::endl;
        return true;
    }
    // ====================================================================
    else if (this->of_conv_ == "both" && potConv && energyConv)
    {
        this->conv_ = true;
        return true;
    }
    else if (this->iter_ >= this->max_iter_)
    {
        return true;
    }
    else
    {
        return false;
    }
}

//
// Print nessecary information
//
void ESolver_OF::printInfo()
{
    if (this->iter_ == 0){
        std::cout << "======================== Running OFDFT ========================" <<  std::endl;
        std::cout << "Iter        Etot(Ha)          Theta      PotNorm     deltaE(Ha)" << std::endl;
        // cout << "======================================== Running OFDFT ========================================" <<  endl;
        // cout << "Iter        Etot(Ha)          Theta       PotNorm        min/max(den)          min/max(dE/dPhi)" << endl;
        // cout << "============================================ OFDFT ========================================" <<  endl;
        // cout << "Iter        Etot(Ha)          Theta       PotNorm        deltaRhoG       deltaRhoR   deltaE" << endl;
    }
    // ============ used to compare with PROFESS3.0 ================
    // double minDen = pelec->charge->rho[0][0];
    // double maxDen = pelec->charge->rho[0][0];
    // double minPot = this->pdEdphi_[0][0];
    // double maxPot = this->pdEdphi_[0][0];
    // for (int i = 0; i < this->pw_rho->nrxx; ++i)
    // {
    //     if (pelec->charge->rho[0][i] < minDen) minDen = pelec->charge->rho[0][i];
    //     if (pelec->charge->rho[0][i] > maxDen) maxDen = pelec->charge->rho[0][i];
    //     if (this->pdEdphi_[0][i] < minPot) minPot = this->pdEdphi_[0][i];
    //     if (this->pdEdphi_[0][i] > maxPot) maxPot = this->pdEdphi_[0][i];
    // }
    std::cout << std::setw(6) << this->iter_
    << std::setw(22) << std::setiosflags(std::ios::scientific) << std::setprecision(12) << this->energy_current_/2.
    << std::setw(12) << std::setprecision(3) << this->theta_[0]
    << std::setw(12) << this->normdLdphi_
    << std::setw(12) << (this->energy_current_ - this->energy_last_)/2. << std::endl;
    // ============ test new convergence criterion =================
    // << setw(12) << this->deltaRhoG
    // << setw(12) << this->deltaRhoR
    // << setw(12) << this->energy_current_ - this->energy_last_ << endl;
    // ============ used to compare with PROFESS3.0 ================
    // << setw(10) << minDen << "/ " << setw(12) << maxDen
    // << setw(10) << minPot << "/ " << setw(10) << maxPot << endl;
    // =============================================================
}

void ESolver_OF::afterOpt(const int istep)
{
    ModuleIO::output_convergence_after_scf(this->conv_, this->pelec->f_en.etot);

    // save charge difference into files for charge extrapolation
    if (GlobalV::CALCULATION != "scf")
    {
        this->CE.save_files(istep,
                            GlobalC::ucell,
#ifdef __MPI
                            this->pw_big,
#endif
                            this->pelec->charge,
                            &this->sf);
    }

    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        if (GlobalV::out_chg == 1)
        {
            std::stringstream ssc;
            ssc << GlobalV::global_out_dir << "SPIN" << is + 1 << "_CHG.cube";
            const double ef_tmp = this->pelec->eferm.get_efval(is);
            ModuleIO::write_rho(
#ifdef __MPI
                pw_big->bz,
                pw_big->nbz,
                pw_rho->nplane,
                pw_rho->startz_current,
#endif
                pelec->charge->rho[is],
                is,
                GlobalV::NSPIN,
                this->iter_,
                ssc.str(),
                pw_rho->nx,
                pw_rho->ny,
                pw_rho->nz,
                ef_tmp,
                &(GlobalC::ucell),
                3);
        }
        
        if (GlobalV::out_pot == 1) // output the effective potential, sunliang 2023-03-16
        {
            int precision = 3; // be consistent with esolver_ks_lcao.cpp
            std::stringstream ssp;
            ssp << GlobalV::global_out_dir << "SPIN" << is + 1 << "_POT.cube";
            ModuleIO::write_potential(
#ifdef __MPI
                this->pw_big->bz,
                this->pw_big->nbz,
                this->pw_rho->nplane,
                this->pw_rho->startz_current,
#endif
                is,
                0,
                ssp.str(),
                this->pw_rho->nx,
                this->pw_rho->ny,
                this->pw_rho->nz,
                this->pelec->pot->get_effective_v(),
                precision);
        }
    }
    if (GlobalV::out_pot == 2) // output the static electronic potential, sunliang 2023-08-11
    {
        int precision = 3;
        std::stringstream ssp;
        ssp << GlobalV::global_out_dir << "/ElecStaticPot.cube";
        ModuleIO::write_elecstat_pot(
#ifdef __MPI
            this->pw_big->bz,
            this->pw_big->nbz,
#endif
            ssp.str(),
            this->pw_rho,
            this->pelec->charge,
            &(GlobalC::ucell),
            this->pelec->pot->get_fixed_v()
        );
    }
}

void ESolver_OF::postprocess()
{

    GlobalV::ofs_running << "\n\n --------------------------------------------" << std::endl;
    GlobalV::ofs_running << std::setprecision(16);
    GlobalV::ofs_running << " !FINAL_ETOT_IS " << this->pelec->f_en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl;
    GlobalV::ofs_running << " --------------------------------------------\n\n" << std::endl;
    // =============== for test ===============
    // if (GlobalV::CAL_FORCE)
    // {
    //     ModuleBase::matrix ff(GlobalC::ucell.nat, 3);
    //     this->cal_Force(ff);
    // }
    // if (GlobalV::CAL_STRESS)
    // {
    //     double unit_transform = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI,3) * 1.0e-8 / 10;
    //     ModuleBase::matrix stress(3,3);
    //     this->cal_Stress(stress);
    //     // stress *= unit_transform;
    //     // cout << "STRESS (GPa)" << endl;
    //     // for (int i = 0; i < 3; ++i)
    //     // {
    //     //     cout << stress(i,0) << "\t"
    //     //         << stress(i, 1) << "\t" << stress(i, 2) << endl;
    //     // }
    // }
}

//
// Get dL/dphi = dL/drho * drho/dphi = (dE/drho - mu) * 2 * ptempPhi and store it in rdLdphi
//
void ESolver_OF::calV(double *ptempPhi, double *rdLdphi)
{
    double **dEdtempPhi = new double*[GlobalV::NSPIN];
    double **tempPhi = new double*[GlobalV::NSPIN];

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        dEdtempPhi[is] = new double[this->pw_rho->nrxx];
        if (is == this->tnSpinFlag_)
        {
            tempPhi[is] = ptempPhi;
        }
        else
        {
            tempPhi[is] = this->pphi_[is];
        }
        for (int ir = 0; ir < this->pw_rho->nrxx; ++ir)
        {
            this->ptempRho_->rho[is][ir] = tempPhi[is][ir] * tempPhi[is][ir];
        }
    }

    if(GlobalV::NSPIN==4) GlobalC::ucell.cal_ux();
    this->pelec->pot->update_from_charge(this->ptempRho_, &GlobalC::ucell);
    ModuleBase::matrix& vr_eff = this->pelec->pot->get_effective_v();

    this->kineticPotential(this->ptempRho_->rho, tempPhi, vr_eff);
    for (int i = 0; i < this->pw_rho->nrxx; ++i)
    {
        dEdtempPhi[this->tnSpinFlag_][i] = vr_eff(this->tnSpinFlag_,i);
    }
    double tempMu = this->cal_mu(ptempPhi, dEdtempPhi[this->tnSpinFlag_], this->nelec_[this->tnSpinFlag_]);
    for (int i = 0; i < this->pw_rho->nrxx; ++i)
    {
        rdLdphi[i] = dEdtempPhi[this->tnSpinFlag_][i] - 2. * tempMu * ptempPhi[i];
    }
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        delete[] dEdtempPhi[is];
    }
    delete[] dEdtempPhi;
    delete[] tempPhi;
}

//
// Calculate dE/dTheta
// dE/dTheta = <dE/dtempPhi|dtempPhi/dTheta>
//           = <dE/dtempPhi|-phi*sin(theta)+d*cos(theta)>
//
void ESolver_OF::caldEdtheta(double **ptempPhi, Charge* tempRho, double *ptheta, double *rdEdtheta)
{
    double *pdPhidTheta = new double[this->pw_rho->nrxx];

    if(GlobalV::NSPIN==4) GlobalC::ucell.cal_ux();
    this->pelec->pot->update_from_charge(tempRho, &GlobalC::ucell);
    ModuleBase::matrix& vr_eff = this->pelec->pot->get_effective_v();

    this->kineticPotential(tempRho->rho, ptempPhi, vr_eff);
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        for (int ir = 0; ir < this->pw_rho->nrxx; ++ir)
        {
            this->pdEdphi_[is][ir] = vr_eff(is,ir);
            pdPhidTheta[ir] = - this->pphi_[is][ir] * sin(ptheta[is]) + this->pdirect_[is][ir] * cos(ptheta[is]);
        }
        rdEdtheta[is] = this->inner_product(this->pdEdphi_[is], pdPhidTheta, this->pw_rho->nrxx, this->dV_);
        Parallel_Reduce::reduce_all(rdEdtheta[is]);
    }
    delete[] pdPhidTheta;
}

//
// Calculate chemical potential mu.
// mu = <dE/dphi|phi> / 2nelec.
//
double ESolver_OF::cal_mu(double *pphi, double *pdEdphi, double nelec)
{
    double mu = this->inner_product(pphi, pdEdphi, this->pw_rho->nrxx, this->dV_);
    Parallel_Reduce::reduce_all(mu);
    mu = mu / (2.0*nelec);
    return mu;
}

// =====================================================================
// NOTE THIS FUNCTION SHOULD BE CALLEDD AFTER POTENTIAL HAS BEEN UPDATED
// =====================================================================
double ESolver_OF::cal_Energy()
{
    this->pelec->cal_energies(2);
    double eKE = this->kineticEnergy(); // kinetic energy
    double ePP = 0.;                    // electron-ion interaction energy
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        ePP += this->inner_product(this->pelec->pot->get_fixed_v(), pelec->charge->rho[is], this->pw_rho->nrxx, this->dV_);
    }
    Parallel_Reduce::reduce_all(ePP);
    this->pelec->f_en.etot += eKE + ePP;
    return this->pelec->f_en.etot;
}

void ESolver_OF::cal_Force(ModuleBase::matrix& force)
{
    Forces<double> ff(GlobalC::ucell.nat);
    ff.cal_force(force, *pelec, this->pw_rho, &this->symm, &sf);
}

void ESolver_OF::cal_Stress(ModuleBase::matrix& stress)
{
    ModuleBase::matrix kinetic_stress_;
    kinetic_stress_.create(3,3);
    this->kinetic_stress(kinetic_stress_);

    OF_Stress_PW ss(this->pelec, this->pw_rho);
    ss.cal_stress(stress, kinetic_stress_, GlobalC::ucell, &this->symm, &sf, &kv);
}
}