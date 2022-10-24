#include "./esolver_of.h"

//-----------temporary-------------------------
#include "../src_pw/global.h"
#include "../module_base/global_function.h"
#include "../module_symmetry/symmetry.h"
// #include "../src_pw/vdwd2.h"
// #include "../src_pw/vdwd3.h"
// #include "../src_pw/vdwd2_parameters.h"
// #include "../src_pw/vdwd3_parameters.h"
#include "../src_pw/pw_complement.h"
#include "../src_pw/structure_factor.h"
#include "../src_pw/symmetry_rho.h"
#include "../src_io/print_info.h"
#include "../src_pw/H_Ewald_pw.h"
#include "../src_pw/electrons.h"
//-----force-------------------
#include "../src_pw/forces.h"
//-----stress------------------
#include "../src_pw/stress_pw.h"
//---------------------------------------------------
#include "module_elecstate/elecstate_pw.h"
#include "module_hamilt/hamilt_pw.h"

namespace ModuleESolver
{

void ESolver_OF::Init(Input &inp, UnitCell_pseudo &ucell)
{
    ESolver_FP::Init(inp, ucell);

    // save necessary parameters
    this->of_kinetic = inp.of_kinetic;
    this->of_method = inp.of_method;
    this->of_conv = inp.of_conv;
    this->of_tole = inp.of_tole;
    this->of_tolp = inp.of_tolp;
    this->maxIter = inp.scf_nmax;

    GlobalC::CHR.cal_nelec();

	if(ucell.atoms[0].xc_func=="HSE"||ucell.atoms[0].xc_func=="PBE0")
	{
        ModuleBase::WARNING_QUIT("esolver_of", "Hybrid functionals are not supported by OFDFT.");
		// XC_Functional::set_xc_type("pbe");
	}
	else
	{
		XC_Functional::set_xc_type(ucell.atoms[0].xc_func);
	}

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SETUP UNITCELL");
    
    // symmetry analysis should be performed every time the cell is changed
    if (ModuleSymmetry::Symmetry::symm_flag)
    {
        GlobalC::symm.analy_sys(ucell, GlobalV::ofs_running);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SYMMETRY");
    }

    // Setup the k points according to symmetry.
    GlobalC::kv.set( GlobalC::symm, GlobalV::global_kpoint_card, GlobalV::NSPIN, ucell.G, ucell.latvec );
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"INIT K-POINTS");

    // print information
    // mohan add 2021-01-30
    Print_Info::setup_parameters(ucell, GlobalC::kv);
    
    // initialize the real-space uniform grid for FFT and parallel
    // distribution of plane waves
    GlobalC::Pgrid.init(GlobalC::rhopw->nx, GlobalC::rhopw->ny, GlobalC::rhopw->nz, GlobalC::rhopw->nplane,
        GlobalC::rhopw->nrxx, GlobalC::bigpw->nbz, GlobalC::bigpw->bz); // mohan add 2010-07-22, update 2011-05-04
    // Calculate Structure factor
    GlobalC::sf.setup_structure_factor(&GlobalC::ucell, GlobalC::rhopw);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT BASIS");

    this->nrxx = this->pw_rho->nrxx;
    this->dV = ucell.omega / GlobalC::rhopw->nxyz; // volume of one point in real space

    //----------------------------------------------------------
    // 1 read in initial data:
    //   a lattice structure:atom_species,atom_positions,lattice vector
    //   b k_points
    //   c pseudopotential
    // 2 setup planeware basis, FFT,structure factor, ...
    // 3 initialize local pseudopotential in G_space
    // 4 initialize charge desity and warefunctios in real space
    //----------------------------------------------------------

    // Inititlize the charge density.
    GlobalC::CHR.allocate(GlobalV::NSPIN, GlobalC::rhopw->nrxx, GlobalC::rhopw->npw);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT CHARGE");

    // Initialize the potential.
    GlobalC::pot.allocate(GlobalC::rhopw->nrxx);

    // Initialize the "wavefunction", which is sqrt(rho)
    this->psi = new psi::Psi<double>(1, GlobalV::NSPIN, this->nrxx);
    this->pphi = new double*[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        this->pphi[is] = this->psi->get_pointer(is);
    }
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT PHI");
 
    //=======================
    // init pseudopotential
    //=======================
    GlobalC::ppcell.init(GlobalC::ucell.ntype);

    //=====================
    // init hamiltonian
    // only allocate in the beginning of ELEC LOOP!
    //=====================
    // GlobalC::hm.hpw.allocate(GlobalC::wf.npwx, GlobalV::NPOL, GlobalC::ppcell.nkb, GlobalC::rhopw->nrxx);

    //=================================
    // initalize local pseudopotential
    //=================================
    GlobalC::ppcell.init_vloc(GlobalC::ppcell.vloc,GlobalC::rhopw);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "LOCAL POTENTIAL");

    //======================================
    // Initalize non local pseudopotential
    //======================================
    GlobalC::ppcell.init_vnl(GlobalC::ucell);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "NON-LOCAL POTENTIAL");

    //=========================================================
    // calculate the total local pseudopotential in real space
    //=========================================================
    GlobalC::pot.init_pot(0, GlobalC::sf.strucFac); //atomic_rho, v_of_rho, set_vrs

    GlobalC::pot.newd();

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT POTENTIAL");

    // Calculate electron numbers
    this->nelec = new double[GlobalV::NSPIN];
    if (GlobalV::NSPIN == 1)
    {
        this->nelec[0] = GlobalC::CHR.nelec;
    }
    else if (GlobalV::NSPIN == 2)
    {
        this->nelec[0] = GlobalC::ucell.magnet.get_nelup();
        this->nelec[1] = GlobalC::ucell.magnet.get_neldw();
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
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT KEDF");
}

void ESolver_OF::Run(int istep, UnitCell_pseudo& ucell)
{
    ModuleBase::timer::tick("ESolver_OF", "Run");
    // get Ewald energy, initial rho and phi if necessary
    this->beforeOpt();
    this->iter = 0;

    while(true)
    {
        // once we get a new rho and phi, update potential
        this->updateV();

        // calculate the energy of new rho and phi
        this->cal_Energy(this->energy_current);

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

    this->afterOpt();

    ModuleBase::timer::tick("ESolver_OF", "Run");
}

// 
// Calculate ewald energy, initialize the rho, phi, theta
// 
void ESolver_OF::beforeOpt()
{
    //calculate ewald energy
    H_Ewald_pw::compute_ewald(GlobalC::ucell, this->pw_rho);

    Symmetry_rho srho;
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        srho.begin(is, GlobalC::CHR, this->pw_rho, GlobalC::Pgrid, GlobalC::symm);
    }

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        if (GlobalC::CHR.init_chg != "file")
        {
            for (int ibs = 0; ibs < this->nrxx; ++ibs)
            {
                // Here we initialize rho to be uniform, 
                // because the rho got by pot.init_pot -> Charge::atomic_rho may contain minus elements.
                GlobalC::CHR.rho[is][ibs] = this->nelec[is]/GlobalC::ucell.omega;
                this->pphi[is][ibs] = sqrt(GlobalC::CHR.rho[is][ibs]);
            }
        }
        else
        {
            for (int ibs = 0; ibs < this->nrxx; ++ibs)
            {
                this->pphi[is][ibs] = sqrt(GlobalC::CHR.rho[is][ibs]);
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
    GlobalC::pot.vr = GlobalC::pot.v_of_rho(GlobalC::CHR.rho, GlobalC::CHR.rho_core); // Hartree + XC
    GlobalC::pot.set_vr_eff(); // Hartree + XC + external
    this->kineticPotential(GlobalC::CHR.rho, this->pphi, GlobalC::pot.vr_eff); // (kinetic + Hartree + XC + external) * 2 * phi
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        for (int ir = 0; ir < this->nrxx; ++ir)
        { 
            this->pdEdphi[is][ir] = GlobalC::pot.vr_eff(is,ir);
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
    double **ptempRho = new double*[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        ptempPhi[is] = new double[this->nrxx];
        ptempRho[is] = new double[this->nrxx];
        for (int ir = 0; ir < this->nrxx; ++ir)
        {
            ptempPhi[is][ir] = this->pphi[is][ir];
            ptempRho[is][ir] = ptempPhi[is][ir] * ptempPhi[is][ir];
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
    this->caldEdtheta(ptempPhi, ptempRho, tempTheta, dEdtheta);
    
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
            this->caldEdtheta(ptempPhi, ptempRho, tempTheta, dEdtheta);
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
//     if (this->iter == 10)
//     {
//         for (int i = -100; i < 100; ++i)
//         {
//             this->theta[0] = 0.001 * i;
//             for (int ir = 0; ir < this->nrxx; ++ir)
//             {
//                 ptempPhi[0][ir] = this->pphi[0][ir] * cos(this->theta[0]) + this->pdirect[0][ir] * sin(this->theta[0]);
//                 ptempRho[0][ir] = ptempPhi[0][ir] * ptempPhi[0][ir];
//             }
//             this->caldEdtheta(ptempPhi, ptempRho, this->theta, dEdtheta);
//             GlobalC::en.calculate_etot();
//             E = GlobalC::en.etot;
//             double eKE = 0.;
//             double ePP = 0.;
//             eKE = this->kineticEnergy();
//             ePP = this->inner_product(GlobalC::pot.vltot, ptempRho[0], this->nrxx, this->dV);
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
            GlobalC::en.calculate_etot();
            E = GlobalC::en.etot;
            eKE = this->kineticEnergy();
            ePP = this->inner_product(GlobalC::pot.vltot, ptempRho[0], this->nrxx, this->dV);
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
                    ptempRho[0][i] = ptempPhi[0][i] * ptempPhi[0][i];
                }
                // get dEdtheta of new tempPhi and tempRho
                this->caldEdtheta(ptempPhi, ptempRho, this->theta, dEdtheta);

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
    //             GlobalC::en.calculate_etot();
    //             E = GlobalC::en.etot;
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
    //                         ptempPhi[is][ir] = this->pphi[is][ir] * cos(tempTheta[is]) + this->pdirect[is][ir] * sin(tempTheta[is]);
    //                         ptempRho[is][ir] = ptempPhi[is][ir] * ptempPhi[is][ir];
    //                     }
    //                 }
    //                 this->caldEdtheta(ptempPhi, ptempRho, tempTheta, dEdtheta);
    //                 dEdalpha = this->inner_product(dEdtheta, pthetaDir, 2, 1);

    //                 if (numDC > 10)
    //                 {
    //                     GlobalV::ofs_warning << "ESolver_OF linesearch: WARNING " << "excedd the max iter number." << endl;
    //                     break;
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
        delete[] ptempRho[is];
    }
    delete[] ptempPhi;
    delete[] ptempRho;
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
        double innerPhiDir = this->inner_product(this->pdirect[0], this->pphi[0], this->nrxx, dV=this->dV);
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
        double normDir = this->inner_product(this->pdirect[0], this->pdirect[0], this->nrxx, dV=this->dV);
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
            double innerPhiDir = this->inner_product(this->pdirect[is], this->pphi[is], this->nrxx, dV=this->dV);
            Parallel_Reduce::reduce_double_all(innerPhiDir);
            for (int i = 0; i < this->nrxx; ++i)
            {
                this->pdirect[is][i] = this->pdirect[is][i] - this->pphi[is][i] * innerPhiDir / this->nelec[is];
            }

            // (2) renormalize direction
            // |d> = |d'> * \sqrt(nelec) / <d'|d'>
            double normDir = this->inner_product(this->pdirect[is], this->pdirect[is], this->nrxx, dV=this->dV);
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
            GlobalC::CHR.rho[is][ir] = this->pphi[is][ir] * this->pphi[is][ir];
        }
    }
    // ============================ for test ===========================
    // if (ModuleSymmetry::Symmetry::symm_flag == 1)
    // {
    //     Symmetry_rho srho;
    //     for (int is = 0; is < GlobalV::NSPIN; is++)
    //     {
    //         srho.begin(is, GlobalC::CHR, this->pw_rho, GlobalC::Pgrid, GlobalC::symm);
    //         for (int ibs = 0; ibs < this->nrxx; ++ibs)
    //         {
    //             this->pphi[is][ibs] = sqrt(GlobalC::CHR.rho[is][ibs]);
    //         }
    //     }
    // }
    // =============test for rho convergence criterion =================
    // for (int is = 0; is < GlobalV::NSPIN; ++is)
    // {
    //     for (int ir = 0; ir < this->nrxx; ++ir)
    //     {
    //         GlobalC::CHR.rho_save[is][ir] = GlobalC::CHR.rho[is][ir];
    //         this->pdeltaRho[is][ir] = GlobalC::CHR.rho[is][ir] - GlobalC::CHR.rho_save[is][ir];
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
        && abs(this->normdLdphi - this->normdLdphi_last) < 1e-10
        && abs(this->normdLdphi - this->normdLdphi_llast) < 1e-10)
        potHold = true;

    if (this->iter >= 3
        && abs(this->energy_current - this->energy_last) < this->of_tole
        && abs(this->energy_current - this->energy_llast) < this->of_tole)
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
    // double minDen = GlobalC::CHR.rho[0][0];
    // double maxDen = GlobalC::CHR.rho[0][0];
    // double minPot = this->pdEdphi[0][0];
    // double maxPot = this->pdEdphi[0][0];
    // for (int i = 0; i < this->nrxx; ++i)
    // {
    //     if (GlobalC::CHR.rho[0][i] < minDen) minDen = GlobalC::CHR.rho[0][i];
    //     if (GlobalC::CHR.rho[0][i] > maxDen) maxDen = GlobalC::CHR.rho[0][i];
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

void ESolver_OF::afterOpt()
{
    if (this->conv)
    {
        GlobalV::ofs_running << "\n charge density convergence is achieved" << std::endl;
        GlobalV::ofs_running << " final etot is " << GlobalC::en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl;
    }
    else
    {
        GlobalV::ofs_running << " convergence has NOT been achieved!" << std::endl;
    }

    if (GlobalC::CHR.out_chg > 0)
    {
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            std::stringstream ssc;
            std::stringstream ss1;
            ssc << GlobalV::global_out_dir << "tmp" << "_SPIN" << is + 1 << "_CHG";
            GlobalC::CHR.write_rho(GlobalC::CHR.rho_save[is], is, iter, ssc.str(), 3);//mohan add 2007-10-17
            ss1 << GlobalV::global_out_dir << "tmp" << "_SPIN" << is + 1 << "_CHG.cube";
            GlobalC::CHR.write_rho_cube(GlobalC::CHR.rho_save[is], is, ss1.str(), 3);
        }
    }
}

void ESolver_OF::postprocess()
{    

    GlobalV::ofs_running << "\n\n --------------------------------------------" << std::endl;
    GlobalV::ofs_running << std::setprecision(16);
    GlobalV::ofs_running << " !FINAL_ETOT_IS " << GlobalC::en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl;
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

    double **tempRho = new double*[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        dEdtempPhi[is] = new double[this->nrxx];
        tempRho[is] = new double[this->nrxx];
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
            tempRho[is][ir] = tempPhi[is][ir] * tempPhi[is][ir];
        }
    }

    GlobalC::pot.vr = GlobalC::pot.v_of_rho(tempRho, GlobalC::CHR.rho_core);
    GlobalC::pot.set_vr_eff();

    this->kineticPotential(tempRho, tempPhi, GlobalC::pot.vr_eff);
    for (int i = 0; i < this->nrxx; ++i)
    {
        dEdtempPhi[this->tnSpinFlag][i] = GlobalC::pot.vr_eff(this->tnSpinFlag,i);
    }
    double tempMu = this->cal_mu(ptempPhi, dEdtempPhi[this->tnSpinFlag], this->nelec[this->tnSpinFlag]);
    for (int i = 0; i < this->nrxx; ++i)
    {
        rdLdphi[i] = dEdtempPhi[this->tnSpinFlag][i] - 2. * tempMu * ptempPhi[i];
    }
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        delete[] dEdtempPhi[is];
        delete[] tempRho[is];
    } 
    delete[] dEdtempPhi;
    delete[] tempRho;
    delete[] tempPhi;
} 

//
// Calculate dE/dTheta
// dE/dTheta = <dE/dtempPhi|dtempPhi/dTheta>
//           = <dE/dtempPhi|-phi*sin(theta)+d*cos(theta)>
//
void ESolver_OF::caldEdtheta(double **ptempPhi, double **ptempRho, double *ptheta, double *rdEdtheta)
{
    double *pdPhidTheta = new double[this->nrxx];

    GlobalC::pot.vr = GlobalC::pot.v_of_rho(ptempRho, GlobalC::CHR.rho_core);
    GlobalC::pot.set_vr_eff();

    this->kineticPotential(ptempRho, ptempPhi, GlobalC::pot.vr_eff);
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        for (int ir = 0; ir < this->nrxx; ++ir)
        {
            this->pdEdphi[is][ir] = GlobalC::pot.vr_eff(is,ir);
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
void ESolver_OF::cal_Energy(double& etot)
{
    GlobalC::en.calculate_etot();
    double eKE = this->kineticEnergy(); // kinetic energy
    double ePP = 0.;                    // electron-ion interaction energy
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        ePP += this->inner_product(GlobalC::pot.vltot, GlobalC::CHR.rho[is], this->nrxx, this->dV);
    }
    Parallel_Reduce::reduce_double_all(ePP);
    GlobalC::en.etot += eKE + ePP;

    this->energy_llast = this->energy_last;
    this->energy_last = this->energy_current;
    this->energy_current = GlobalC::en.etot;
}


void ESolver_OF::cal_Force(ModuleBase::matrix& force)
{
    Forces ff;
    ff.init(force);
}

void ESolver_OF::cal_Stress(ModuleBase::matrix& stress)
{
    Stress_PW ss;
    ss.cal_stress(stress);
    if (this->of_kinetic == "tf")
    {    
        this->tf.get_stress(GlobalC::ucell.omega);
        stress += this->tf.stress;
    }
    else if (this->of_kinetic == "vw")
    {
        this->vw.get_stress(this->pphi, this->pw_rho);
        stress += this->vw.stress;
    }
    else if (this->of_kinetic == "wt")
    {
        this->tf.get_stress(GlobalC::ucell.omega);
        this->vw.get_stress(this->pphi, this->pw_rho);
        this->wt.get_stress(GlobalC::ucell.omega, GlobalC::CHR.rho, this->pw_rho, GlobalV::of_vw_weight);
        stress += this->tf.stress + this->vw.stress + this->wt.stress;
    }
    else if (this->of_kinetic == "tf+")
    {
        this->tf.get_stress(GlobalC::ucell.omega);
        this->vw.get_stress(this->pphi, this->pw_rho);
        stress += this->tf.stress + this->vw.stress;
    }
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
    return kinetic;
}
}