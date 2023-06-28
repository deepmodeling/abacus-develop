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
    this->of_kinetic = inp.of_kinetic;
    this->of_method = inp.of_method;
    this->of_conv = inp.of_conv;
    this->of_tole = inp.of_tole;
    this->of_tolp = inp.of_tolp;
    this->maxIter = inp.scf_nmax;

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

    this->nrxx = this->pw_rho->nrxx;
    this->dV = ucell.omega / this->pw_rho->nxyz; // volume of one point in real space

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
    this->psi = new psi::Psi<double>(1, GlobalV::NSPIN, this->nrxx);
    ModuleBase::Memory::record("OFDFT::Psi", sizeof(double) * GlobalV::NSPIN * this->nrxx);
    this->pphi = new double*[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        this->pphi[is] = this->psi->get_pointer(is);
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
    GlobalC::ppcell.init_vnl(GlobalC::ucell);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "NON-LOCAL POTENTIAL");

    GlobalC::ppcell.cal_effective_D();

    if(this->pelec == nullptr)
    {
        this->pelec = new elecstate::ElecState((Charge*)(&chr), this->pw_rho, pw_big);
    }

    this->pelec->charge->allocate(GlobalV::NSPIN);
    this->pelec->omega = GlobalC::ucell.omega;

    this->pelec->pot = new elecstate::Potential(pw_rho,
                                                &GlobalC::ucell,
                                                &GlobalC::ppcell.vloc,
                                                &sf,
                                                &(this->pelec->f_en.etxc),
                                                &(this->pelec->f_en.vtxc));
    //There is no Operator in ESolver_OF, register Potentials here!
    std::vector<string> pot_register_in;
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

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT POTENTIAL");

    // Calculate electron numbers
    this->nelec = new double[GlobalV::NSPIN];
    if (GlobalV::NSPIN == 1)
    {
        this->nelec[0] = GlobalV::nelec;
    }
    else if (GlobalV::NSPIN == 2)
    {
        //in fact, nelec_spin will not be used anymore
        this->pelec->init_nelec_spin();
        this->nelec[0] = this->pelec->nelec_spin[0];
        this->nelec[1] = this->pelec->nelec_spin[1];
    }

    // ================================
    // Initialize optimization methods
    // ================================
    if (this->of_method == "tn")
    {
        this->opt_tn.allocate(this->nrxx);
        this->opt_tn.setPara(this->dV);
    }
    else if (this->of_method == "cg1" || this->of_method == "cg2")
    {
        this->opt_cg.allocate(this->nrxx);
        this->opt_cg.setPara(this->dV);
        this->opt_dcsrch.set_paras(1e-4,1e-2);
    }
    else if (this->of_method == "bfgs")
    {
        ModuleBase::WARNING_QUIT("esolver_of", "BFGS is not supported now.");
        return;
    }

    // optimize theta if nspin=2
    if (GlobalV::NSPIN == 2)
    {
        this->opt_cg_mag = new ModuleBase::Opt_CG;
        this->opt_cg_mag->allocate(GlobalV::NSPIN);
    }

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT OPTIMIZATION");

    // =============================================
    // Initalize chemical potential, step length, ...
    // =============================================
    this->mu = new double[GlobalV::NSPIN];
    this->theta = new double[GlobalV::NSPIN];
    this->pdLdphi = new double*[GlobalV::NSPIN];
    this->pdEdphi = new double*[GlobalV::NSPIN];
    this->pdirect = new double*[GlobalV::NSPIN];
    this->precipDir = new std::complex<double> *[GlobalV::NSPIN];
    // this->pdeltaRhoHar = new double[this->pw_rho->nrxx];

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        this->pdLdphi[is] = new double[this->nrxx];
        this->pdEdphi[is] = new double[this->nrxx];
        this->pdirect[is] = new double[this->nrxx];
        this->precipDir[is] = new std::complex<double>[pw_rho->npw];
    }

    // ===================================
    // Initialize KEDF
    // ===================================
    this->tf.set_para(this->nrxx, this->dV, GlobalV::of_tf_weight);
    this->vw.set_para(this->nrxx, this->dV, GlobalV::of_vw_weight);
    this->wt.set_para(this->nrxx, this->dV, GlobalV::of_wt_alpha, GlobalV::of_wt_beta, this->nelec[0], GlobalV::of_tf_weight, GlobalV::of_vw_weight, GlobalV::of_read_kernel, GlobalV::of_kernel_file, this->pw_rho);
    this->lkt.set_para(this->nrxx, this->dV, GlobalV::of_lkt_a);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT KEDF");

    // Initialize charge extrapolation
    CE.Init_CE(GlobalC::ucell.nat);
    delete this->ptempRho;
    this->ptempRho = new Charge();
    this->ptempRho->set_rhopw(this->pw_rho);
    this->ptempRho->allocate(GlobalV::NSPIN);
}

void ESolver_OF::init_after_vc(Input &inp, UnitCell &ucell)
{
    ModuleBase::timer::tick("ESolver_OF", "init_after_vc");

    ESolver_FP::init_after_vc(inp,ucell);

    GlobalC::ppcell.init_vnl(GlobalC::ucell);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"NON-LOCAL POTENTIAL");
}

void ESolver_OF::Run(int istep, UnitCell& ucell)
{
    ModuleBase::timer::tick("ESolver_OF", "Run");
    // get Ewald energy, initial rho and phi if necessary
    this->beforeOpt(istep);
    this->iter = 0;

    while(true)
    {
        // once we get a new rho and phi, update potential
        this->updateV();

        // calculate the energy of new rho and phi
        this->energy_llast = this->energy_last;
        this->energy_last = this->energy_current;
        this->energy_current = this->cal_Energy();

        // print neccesary information
        this->printInfo();

        // check if the job is done
        if (this->checkExit()) break;

        // find the optimization direction and step lenghth theta according to the potential
        this->solveV();

        // update the rho and phi based on the direction and theta
        this->updateRho();

        this->iter++;
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
            for (int ibs = 0; ibs < this->nrxx; ++ibs)
            {
                // Here we initialize rho to be uniform,
                // because the rho got by pot.init_pot -> Charge::atomic_rho may contain minus elements.
                pelec->charge->rho[is][ibs] = this->nelec[is]/GlobalC::ucell.omega;
                this->pphi[is][ibs] = sqrt(pelec->charge->rho[is][ibs]);
            }
        }
        else
        {
            for (int ibs = 0; ibs < this->nrxx; ++ibs)
            {
                this->pphi[is][ibs] = sqrt(pelec->charge->rho[is][ibs]);
            }
        }
    }

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        this->mu[is] = 0;
        this->theta[is] = 0.;
        ModuleBase::GlobalFunc::ZEROS(this->pdLdphi[is], this->nrxx);
        ModuleBase::GlobalFunc::ZEROS(this->pdEdphi[is], this->nrxx);
        ModuleBase::GlobalFunc::ZEROS(this->pdirect[is], this->nrxx);
    }
    if (GlobalV::NSPIN == 1)
    {
        this->theta[0] = 0.2;
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
    this->kineticPotential(pelec->charge->rho, this->pphi, this->pelec->pot->get_effective_v()); // (kinetic + Hartree + XC + external) * 2 * phi
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        const double* vr_eff = this->pelec->pot->get_effective_v(is);
        for (int ir = 0; ir < this->nrxx; ++ir)
        {
            this->pdEdphi[is][ir] = vr_eff[ir];
        }
        this->mu[is] = this->cal_mu(this->pphi[is], this->pdEdphi[is], this->nelec[is]);

        for (int ir = 0; ir < this->nrxx; ++ir)
        {
            this->pdLdphi[is][ir] = this->pdEdphi[is][ir] - 2. * this->mu[is] * this->pphi[is][ir];
        }
    }

    // (2) get the norm of dLdphi
    // ===== temporary solution of potential convergence when of_full_pw = 0 =====
    this->normdLdphi_llast = this->normdLdphi_last;
    this->normdLdphi_last = this->normdLdphi;
    // ===========================================================================
    this->normdLdphi = 0.;

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
       this->normdLdphi += this->inner_product(this->pdLdphi[is], this->pdLdphi[is], this->nrxx, 1);
    }
    Parallel_Reduce::reduce_double_all(this->normdLdphi);
    this->normdLdphi = sqrt(this->normdLdphi/this->pw_rho->nxyz/GlobalV::NSPIN);
}

//
// Get optimization direction d and step theta
//
void ESolver_OF::solveV()
{
    // (1) get |d0> with optimization algorithm
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        if (this->of_method == "tn")
        {
            this->tnSpinFlag = is;
            opt_tn.next_direct(this->pphi[is], this->pdLdphi[is], flag, this->pdirect[is], this, &ESolver_OF::calV);
        }
        else if (this->of_method == "cg1")
        {
            opt_cg.next_direct(this->pdLdphi[is], 1, this->pdirect[is]);
        }
        else if (this->of_method == "cg2")
        {
            opt_cg.next_direct(this->pdLdphi[is], 2, this->pdirect[is]);
        }
        else if (this->of_method == "bfgs")
        {
            return;
        }
        else
        {
            ModuleBase::WARNING_QUIT("ESolver_OF", "of_method must be one of CG, TN, or BFGS.");
        }
    }
    // initialize tempPhi and tempRho used in line search
    double **ptempPhi = new double*[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        ptempPhi[is] = new double[this->nrxx];
        for (int ir = 0; ir < this->nrxx; ++ir)
        {
            ptempPhi[is][ir] = this->pphi[is][ir];
            this->ptempRho->rho[is][ir] = ptempPhi[is][ir] * ptempPhi[is][ir];
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
    this->caldEdtheta(ptempPhi, this->ptempRho, tempTheta, dEdtheta);

    // Assert dEdtheta(theta = 0) < 0, otherwise line search will not work.
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        if (dEdtheta[is] > dEdthetaThre)
        {
            cout << "dEdtheta    " << dEdtheta[is] << endl;
            ModuleBase::WARNING_QUIT("esolver_of.cpp", "dE/dtheta is too large.");
        }
        else if (dEdtheta[is] > 0)
        {
            GlobalV::ofs_warning << "ESolver_OF: WARNING " << "dEdphi > 0, replace direct with steepest descent method." << endl;
            for (int ir = 0; ir < this->nrxx; ++ir)
            {
                this->pdirect[is][ir] = - this->pdLdphi[is][ir];
            }
            this->getNextDirect();
            this->caldEdtheta(ptempPhi, this->ptempRho, tempTheta, dEdtheta);
            if (dEdtheta[is] > dEdthetaThre)
            {
                cout << "dEdtheta    " << dEdtheta[is] << endl;
                ModuleBase::WARNING_QUIT("esolver_of.cpp", "dE/dtheta is too large.");
            }
            else if (dEdtheta[is] > 0)
            {
                GlobalV::ofs_warning << "ESolver_OF: WARNING " << "when use steepest dencent method, dEdphi > 0, so we might get minimum." << endl;
            }
        }
    }
    delete[] tempTheta;

    // // ======================== for test ============================
    //     if (this->iter == 0)
    //     {
    //         for (int i = -100; i < 100; ++i)
    //         {
    //             this->theta[0] = 0.001 * i;
    //             for (int ir = 0; ir < this->nrxx; ++ir)
    //             {
    //                 ptempPhi[0][ir] = this->pphi[0][ir] * cos(this->theta[0]) + this->pdirect[0][ir] *
    //                 sin(this->theta[0]); ptempRho->rho[0][ir] = ptempPhi[0][ir] * ptempPhi[0][ir];
    //             }
    //             this->caldEdtheta(ptempPhi, ptempRho, this->theta, dEdtheta);
    //             this->pelec->f_en.calculate_etot(this->pw_rho->nrxx, this->pw_rho->nxyz);
    //             E = this->pelec->f_en.etot;
    //             double eKE = 0.;
    //             double ePP = 0.;
    //             eKE = this->kineticEnergy();
    //             ePP = this->inner_product(this->pelec->pot->get_fixed_v(), ptempRho->rho[0], this->nrxx, this->dV);
    //             // ePP = this->inner_product(GlobalC::pot.vltot, ptempRho[0], this->nrxx, this->dV);
    //             Parallel_Reduce::reduce_double_all(ePP);
    //             E += eKE + ePP;
    //             GlobalV::ofs_warning << i << "    " << dEdtheta[0] << "    " << E << endl;
    //             if (this->theta[0] == 0) cout << "dEdtheta    " << dEdtheta[0]<< endl;
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
        strcpy(this->task, "START");
        while (true)
        {
            // update energy
            this->pelec->cal_energies(2);
            E = this->pelec->f_en.etot;
            eKE = this->kineticEnergy();
            ePP = this->inner_product(this->pelec->pot->get_fixed_v(), this->ptempRho->rho[0], this->nrxx, this->dV);
            Parallel_Reduce::reduce_double_all(ePP);
            E += eKE + ePP;

            // line search to update theta[0]
            this->opt_dcsrch.dcSrch(E, dEdtheta[0], this->theta[0], this->task);
            numDC++;

            // decide what to do next according to the output of line search
            if (strncmp(this->task, "FG", 2) == 0) // continue line search
            {
                // update tempPhi and tempRho
                for (int i = 0; i < this->nrxx; ++i)
                {
                    ptempPhi[0][i] = this->pphi[0][i] * cos(this->theta[0]) + this->pdirect[0][i] * sin(this->theta[0]);
                    this->ptempRho->rho[0][i] = ptempPhi[0][i] * ptempPhi[0][i];
                }
                // get dEdtheta of new tempPhi and tempRho
                this->caldEdtheta(ptempPhi, this->ptempRho, this->theta, dEdtheta);

                if (numDC > this->maxDCsrch)
                {
                    GlobalV::ofs_warning << "ESolver_OF linesearch: WARNING " << "excedd the max iter number." << endl;
                    break;
                }
            }
            else if (strncmp(this->task, "CO", 2) == 0) // convergence achieved
            {
                break;
            }
            else if (strncmp(this->task, "WA", 2) == 0) // warning of line search
            {
                GlobalV::ofs_warning << "ESolver_OF linesearch: WARNING " << this->task << std::endl;
                cout << this->task << endl;
                break;
            }
            else if (strncmp(this->task, "ER", 2) == 0) // ERROR in line search
            {
                GlobalV::ofs_warning << "ESolver_OF linesearch: ERROR " << this->task << std::endl;
                cout << this->task << endl;
                break;
            }
        }
    }
    else if (GlobalV::NSPIN == 2)
    {
        ModuleBase::WARNING_QUIT("esolver_of", "Sorry, SPIN2 case is not supported by OFDFT for now.");
    // ========================== Under testing ==========================
    //     this->opt_cg_mag->refresh();

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
    //         this->opt_cg_mag->next_direct(dEdtheta, 1, pthetaDir);

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
        //         this->opt_dcsrch.set_paras(1e-4, 1e-2, 1e-12, 0., ModuleBase::PI/maxThetaDir);
        //         strcpy(this->task, "START");
        //         numDC = 0;
        //         while(true)
        //         {
        //             this->pelec->f_en.calculate_etot(this->pw_rho->nrxx, this->pw_rho->nxyz);
        //             E = this->pelec->f_en.etot;
        //             eKE = this->kineticEnergy();
        //             ePP = 0.;
        //             for (int is = 0; is < GlobalV::NSPIN; ++is) {
        //                 ePP += this->inner_product(GlobalC::pot.vltot, ptempRho[is], this->nrxx, this->dV);
        //             }
        //             Parallel_Reduce::reduce_double_all(ePP);
        //             E += eKE + ePP;
        //             this->opt_dcsrch.dcSrch(E, dEdalpha, thetaAlpha, this->task);
        //             numDC++;

        //             if (strncmp(this->task, "FG", 2) == 0)
        //             {
        //                 for (int is = 0; is < GlobalV::NSPIN; ++is)
        //                 {
        //                     tempTheta[is] = this->theta[is] + thetaAlpha * pthetaDir[is];
        //                     for (int ir = 0; ir < this->nrxx; ++ir)
        //                     {
        //                         ptempPhi[is][ir] = this->pphi[is][ir] * cos(tempTheta[is]) + this->pdirect[is][ir] *
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
        //             else if (strncmp(this->task, "CO", 2) == 0)
        //             {
        //                 break;
        //             }
        //             else if (strncmp(this->task, "WA", 2) == 0)
        //             {
        //                 GlobalV::ofs_warning << "ESolver_OF linesearch: WARNING " << this->task << std::endl;
        //                 cout << this->task << endl;
        //                 break;
        //             }
        //             else if (strncmp(this->task, "ER", 2) == 0)
        //             {
        //                 GlobalV::ofs_warning << "ESolver_OF linesearch: ERROR " << this->task << std::endl;
        //                 cout << this->task << endl;
        //                 break;
        //             }
        //         }

        //         for (int is = 0; is < GlobalV::NSPIN; ++is) this->theta[is] += thetaAlpha * pthetaDir[is];
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
            pw_rho->real2recip(this->pdirect[is], this->precipDir[is]);
            pw_rho->recip2real(this->precipDir[is], this->pdirect[is]);
        }
    }

    if (GlobalV::NSPIN == 1)
    {
        double tempTheta = 0; // tempTheta = |d'|/|d0 + phi|, theta = min(theta, tempTheta)

        // (1) make direction orthogonal to phi
        // |d'> = |d0> - |phi><phi|d0>/nelec
        double innerPhiDir = this->inner_product(this->pdirect[0], this->pphi[0], this->nrxx, this->dV);
        Parallel_Reduce::reduce_double_all(innerPhiDir);
        for (int i = 0; i < this->nrxx; ++i)
        {
            tempTheta += pow(this->pdirect[0][i] + this->pphi[0][i], 2);
            this->pdirect[0][i] = this->pdirect[0][i] - this->pphi[0][i] * innerPhiDir / this->nelec[0];
        }
        Parallel_Reduce::reduce_double_all(tempTheta);
        tempTheta = sqrt(tempTheta);

        // (2) renormalize direction
        // |d> = |d'> * \sqrt(nelec) / <d'|d'>
        double normDir = this->inner_product(this->pdirect[0], this->pdirect[0], this->nrxx, this->dV);
        Parallel_Reduce::reduce_double_all(normDir);
        normDir = sqrt(normDir);
        for (int i = 0; i < this->nrxx; ++i)
        {
            this->pdirect[0][i] = sqrt(this->nelec[0]) * this->pdirect[0][i] / normDir;
        }

        tempTheta = normDir/tempTheta;
        this->theta[0] = min(this->theta[0], tempTheta);
    }
    else if (GlobalV::NSPIN == 2) // theta = 0
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            // (1) make direction orthogonal to phi
            // |d'> = |d0> - |phi><phi|d0>/nelec
            double innerPhiDir = this->inner_product(this->pdirect[is], this->pphi[is], this->nrxx, this->dV);
            Parallel_Reduce::reduce_double_all(innerPhiDir);
            for (int i = 0; i < this->nrxx; ++i)
            {
                this->pdirect[is][i] = this->pdirect[is][i] - this->pphi[is][i] * innerPhiDir / this->nelec[is];
            }

            // (2) renormalize direction
            // |d> = |d'> * \sqrt(nelec) / <d'|d'>
            double normDir = this->inner_product(this->pdirect[is], this->pdirect[is], this->nrxx, this->dV);
            Parallel_Reduce::reduce_double_all(normDir);
            normDir = sqrt(normDir);
            for (int i = 0; i < this->nrxx; ++i)
            {
                this->pdirect[is][i] = sqrt(this->nelec[is]) * this->pdirect[is][i] / normDir;
            }
            this->theta[is] = 0.;
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
        for (int ir = 0; ir < this->nrxx; ++ir)
        {
            this->pphi[is][ir] = this->pphi[is][ir] * cos(this->theta[is]) + this->pdirect[is][ir] * sin(this->theta[is]);
            pelec->charge->rho[is][ir] = this->pphi[is][ir] * this->pphi[is][ir];
        }
    }
    // ============================ for test ===========================
    // if (ModuleSymmetry::Symmetry::symm_flag == 1)
    // {
    //     Symmetry_rho srho;
    //     for (int is = 0; is < GlobalV::NSPIN; is++)
    //     {
    //         srho.begin(is, pelec->charge, this->pw_rho, GlobalC::Pgrid, this->symm);
    //         for (int ibs = 0; ibs < this->nrxx; ++ibs)
    //         {
    //             this->pphi[is][ibs] = sqrt(pelec->charge->rho[is][ibs]);
    //         }
    //     }
    // }
    // =============test for rho convergence criterion =================
    // for (int is = 0; is < GlobalV::NSPIN; ++is)
    // {
    //     for (int ir = 0; ir < this->nrxx; ++ir)
    //     {
    //         pelec->charge->rho_save[is][ir] = pelec->charge->rho[is][ir];
    //         this->pdeltaRho[is][ir] = pelec->charge->rho[is][ir] - pelec->charge->rho_save[is][ir];
    //         this->deltaRhoR += abs(this->pdeltaRho[is][ir]);

    //     }
    //     this->pw_rho->real2recip(this->pdeltaRho[is], this->precipDir[is]);
    //     for (int ig = 0; ig < this->pw_rho->npw; ++ig)
    //     {
    //         if (this->pw_rho->gg[ig] != 0.)
    //             this->precipDir[is][ig] = this->precipDir[is][ig] / this->pw_rho->gg[ig] / this->pw_rho->tpiba2 * 4. * M_PI;
    //         else
    //             this->precipDir[is][ig] = 0.;
    //     }
    //     this->pw_rho->recip2real(this->precipDir[is], this->pdeltaRhoHar);
    //     this->deltaRhoG = this->inner_product(this->pdeltaRho[is], this->pdeltaRhoHar, this->nrxx, this->dV);
    // }
    // Parallel_Reduce::reduce_double_all(this->deltaRhoR);
    // Parallel_Reduce::reduce_double_all(this->deltaRhoG);
    // this->deltaRhoR *= this->dV;
    // this->deltaRhoG /= 2.;
}

//
// Check convergence, return ture if converge or iter >= maxIter.
//
bool ESolver_OF::checkExit()
{
    this->conv = false;
    bool potConv = false;
    bool potHold = false; // if normdLdphi nearly remains unchanged
    bool energyConv = false;

    if (this->normdLdphi < this->of_tolp)
        potConv = true;
    if (this->iter >= 3
        && std::abs(this->normdLdphi - this->normdLdphi_last) < 1e-10
        && std::abs(this->normdLdphi - this->normdLdphi_llast) < 1e-10)
        potHold = true;

    if (this->iter >= 3
        && std::abs(this->energy_current - this->energy_last) < this->of_tole
        && std::abs(this->energy_current - this->energy_llast) < this->of_tole)
        energyConv = true;

    if (this->of_conv == "energy" && energyConv)
    {
        this->conv = true;
        return true;
    }
    else if (this->of_conv == "potential" && potConv)
    {
        this->conv = true;
        return true;
    // ============ temporary solution of potential convergence ===========
    }
    else if (this->of_conv == "potential" && potHold)
    {
        GlobalV::ofs_warning << "ESolver_OF WARNING: " <<
        "The convergence of potential has not been reached, but the norm of potential nearly remains unchanged, set of_full_pw = 1 may work." << endl;
        return true;
    }
    // ====================================================================
    else if (this->of_conv == "both" && potConv && energyConv)
    {
        this->conv = true;
        return true;
    }
    else if (this->iter >= this->maxIter)
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
    if (this->iter == 0){
        cout << "======================== Running OFDFT ========================" <<  endl;
        cout << "Iter        Etot(Ha)          Theta      PotNorm     deltaE(Ha)" << endl;
        // cout << "======================================== Running OFDFT ========================================" <<  endl;
        // cout << "Iter        Etot(Ha)          Theta       PotNorm        min/max(den)          min/max(dE/dPhi)" << endl;
        // cout << "============================================ OFDFT ========================================" <<  endl;
        // cout << "Iter        Etot(Ha)          Theta       PotNorm        deltaRhoG       deltaRhoR   deltaE" << endl;
    }
    // ============ used to compare with PROFESS3.0 ================
    // double minDen = pelec->charge->rho[0][0];
    // double maxDen = pelec->charge->rho[0][0];
    // double minPot = this->pdEdphi[0][0];
    // double maxPot = this->pdEdphi[0][0];
    // for (int i = 0; i < this->nrxx; ++i)
    // {
    //     if (pelec->charge->rho[0][i] < minDen) minDen = pelec->charge->rho[0][i];
    //     if (pelec->charge->rho[0][i] > maxDen) maxDen = pelec->charge->rho[0][i];
    //     if (this->pdEdphi[0][i] < minPot) minPot = this->pdEdphi[0][i];
    //     if (this->pdEdphi[0][i] > maxPot) maxPot = this->pdEdphi[0][i];
    // }
    cout << setw(6) << this->iter
    << setw(22) << setiosflags(ios::scientific) << setprecision(12) << this->energy_current/2.
    << setw(12) << setprecision(3) << this->theta[0]
    << setw(12) << this->normdLdphi
    << setw(12) << (this->energy_current - this->energy_last)/2. << endl;
    // ============ test new convergence criterion =================
    // << setw(12) << this->deltaRhoG
    // << setw(12) << this->deltaRhoR
    // << setw(12) << this->energy_current - this->energy_last << endl;
    // ============ used to compare with PROFESS3.0 ================
    // << setw(10) << minDen << "/ " << setw(12) << maxDen
    // << setw(10) << minPot << "/ " << setw(10) << maxPot << endl;
    // =============================================================
}

void ESolver_OF::afterOpt(const int istep)
{
    ModuleIO::output_convergence_after_scf(this->conv, this->pelec->f_en.etot);

    // save charge difference into files for charge extrapolation
    this->CE.save_files(istep,
                        GlobalC::ucell,
#ifdef __MPI
                        this->pw_big,
#endif
                        this->pelec->charge,
                        &this->sf);

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
                iter,
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
                pw_big->bz,
                pw_big->nbz,
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
        dEdtempPhi[is] = new double[this->nrxx];
        if (is == this->tnSpinFlag)
        {
            tempPhi[is] = ptempPhi;
        }
        else
        {
            tempPhi[is] = this->pphi[is];
        }
        for (int ir = 0; ir < this->nrxx; ++ir)
        {
            this->ptempRho->rho[is][ir] = tempPhi[is][ir] * tempPhi[is][ir];
        }
    }

    if(GlobalV::NSPIN==4) GlobalC::ucell.cal_ux();
    this->pelec->pot->update_from_charge(this->ptempRho, &GlobalC::ucell);
    ModuleBase::matrix& vr_eff = this->pelec->pot->get_effective_v();

    this->kineticPotential(this->ptempRho->rho, tempPhi, vr_eff);
    for (int i = 0; i < this->nrxx; ++i)
    {
        dEdtempPhi[this->tnSpinFlag][i] = vr_eff(this->tnSpinFlag,i);
    }
    double tempMu = this->cal_mu(ptempPhi, dEdtempPhi[this->tnSpinFlag], this->nelec[this->tnSpinFlag]);
    for (int i = 0; i < this->nrxx; ++i)
    {
        rdLdphi[i] = dEdtempPhi[this->tnSpinFlag][i] - 2. * tempMu * ptempPhi[i];
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
    double *pdPhidTheta = new double[this->nrxx];

    if(GlobalV::NSPIN==4) GlobalC::ucell.cal_ux();
    this->pelec->pot->update_from_charge(tempRho, &GlobalC::ucell);
    ModuleBase::matrix& vr_eff = this->pelec->pot->get_effective_v();

    this->kineticPotential(tempRho->rho, ptempPhi, vr_eff);
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        for (int ir = 0; ir < this->nrxx; ++ir)
        {
            this->pdEdphi[is][ir] = vr_eff(is,ir);
            pdPhidTheta[ir] = - this->pphi[is][ir] * sin(ptheta[is]) + this->pdirect[is][ir] * cos(ptheta[is]);
        }
        rdEdtheta[is] = this->inner_product(this->pdEdphi[is], pdPhidTheta, this->nrxx, this->dV);
        Parallel_Reduce::reduce_double_all(rdEdtheta[is]);
    }
    delete[] pdPhidTheta;
}

//
// Calculate chemical potential mu.
// mu = <dE/dphi|phi> / 2nelec.
//
double ESolver_OF::cal_mu(double *pphi, double *pdEdphi, double nelec)
{
    double mu = this->inner_product(pphi, pdEdphi, this->nrxx, this->dV);
    Parallel_Reduce::reduce_double_all(mu);
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
        ePP += this->inner_product(this->pelec->pot->get_fixed_v(), pelec->charge->rho[is], this->nrxx, this->dV);
    }
    Parallel_Reduce::reduce_double_all(ePP);
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
    ModuleBase::matrix kinetic_stress;
    kinetic_stress.create(3,3);
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            kinetic_stress(i,j) = 0.0;
        }
    }

    if (this->of_kinetic == "tf")
    {
        this->tf.get_stress(GlobalC::ucell.omega);
        kinetic_stress += this->tf.stress;
    }
    else if (this->of_kinetic == "vw")
    {
        this->vw.get_stress(this->pphi, this->pw_rho);
        kinetic_stress += this->vw.stress;
    }
    else if (this->of_kinetic == "wt")
    {
        this->tf.get_stress(GlobalC::ucell.omega);
        this->vw.get_stress(this->pphi, this->pw_rho);
        this->wt.get_stress(GlobalC::ucell.omega, pelec->charge->rho, this->pw_rho, GlobalV::of_vw_weight);
        kinetic_stress += this->tf.stress + this->vw.stress + this->wt.stress;
    }
    else if (this->of_kinetic == "tf+")
    {
        this->tf.get_stress(GlobalC::ucell.omega);
        this->vw.get_stress(this->pphi, this->pw_rho);
        kinetic_stress += this->tf.stress + this->vw.stress;
    }
    else if (this->of_kinetic == "lkt")
    {
        this->lkt.get_stress(GlobalC::ucell.omega, pelec->charge->rho, this->pw_rho);
        this->vw.get_stress(pelec->charge->rho, this->pw_rho);
        kinetic_stress += this->lkt.stress + this->vw.stress;
    }

    OF_Stress_PW ss(this->pelec, this->pw_rho);
    ss.cal_stress(stress, kinetic_stress, GlobalC::ucell, &this->symm, &sf, &kv);
}

// Calculated kinetic potential and plus it to &rpot, return (rpot + kietic potential) * 2 * pphiInpt
void ESolver_OF::kineticPotential(double **prho, double **pphiInpt, ModuleBase::matrix &rpot)
{
    if (this->of_kinetic == "tf")
    {
        this->tf.tf_potential(prho, rpot);
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            for (int ir = 0; ir < this->nrxx; ++ir)
            {
                rpot(is,ir) *= 2.0 * pphiInpt[is][ir];
            }
        }
    }
    else if (this->of_kinetic == "vw")
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            for (int ir = 0; ir < this->nrxx; ++ir)
            {
                rpot(is,ir) *= 2.0 * pphiInpt[is][ir];
            }
        }
        this->vw.vW_potential(pphiInpt, this->pw_rho, rpot);
    }
    else if (this->of_kinetic == "wt")
    {
        this->tf.tf_potential(prho, rpot);
        this->wt.WT_potential(prho, this->pw_rho, rpot);
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            for (int ir = 0; ir < this->nrxx; ++ir)
            {
                rpot(is,ir) *= 2.0 * pphiInpt[is][ir];
            }
        }
        this->vw.vW_potential(pphiInpt, this->pw_rho, rpot);
    }
    else if (this->of_kinetic == "tf+")
    {
        this->tf.tf_potential(prho, rpot);
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            for (int ir = 0; ir < this->nrxx; ++ir)
            {
                rpot(is,ir) *= 2.0 * pphiInpt[is][ir];
            }
        }
        this->vw.vW_potential(pphiInpt, this->pw_rho, rpot);
    }
    else if (this->of_kinetic == "lkt")
    {
        this->lkt.lkt_potential(prho, this->pw_rho, rpot);
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            for (int ir = 0; ir < this->nrxx; ++ir)
            {
                rpot(is,ir) *= 2.0 * pphiInpt[is][ir];
            }
        }
        this->vw.vW_potential(pphiInpt, this->pw_rho, rpot);
    }
}

// Return the kinetic energy
double ESolver_OF::kineticEnergy()
{
    double kinetic = 0.;
    if (this->of_kinetic == "tf")
    {
        kinetic += this->tf.TFenergy;
    }
    else if (this->of_kinetic == "vw")
    {
        kinetic += this->vw.vWenergy;
    }
    else if (this->of_kinetic == "wt")
    {
        kinetic += this->tf.TFenergy + this->vw.vWenergy + this->wt.WTenergy;
    }
    else if (this->of_kinetic == "tf+")
    {
        kinetic += this->tf.TFenergy + this->vw.vWenergy;
    }
    else if (this->of_kinetic == "lkt")
    {
        kinetic += this->lkt.LKTenergy + this->vw.vWenergy;
    }
    return kinetic;
}
}