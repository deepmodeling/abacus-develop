#include "esolver_of.h"
#include "module_base/memory.h"
#include "module_base/formatter.h"
#include "module_elecstate/potentials/efield.h"
#include "module_elecstate/potentials/gatefield.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace ModuleESolver
{

/**
 * @brief Initialize this->pelec, as well as this->pelec->pot
 * 
 * @param ucell 
 */
void ESolver_OF::init_elecstate(UnitCell &ucell)
{
    delete this->pelec;
    this->pelec = new elecstate::ElecState((Charge*)(&chr), this->pw_rho, pw_big);

    this->pelec->charge->allocate(GlobalV::NSPIN);
    this->pelec->omega = ucell.omega;

    delete this->pelec->pot;
    this->pelec->pot = new elecstate::Potential(this->pw_rhod,
                                                this->pw_rho,
                                                &ucell,
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
}

/**
 * @brief Allocate the arrays, as well as this->psi_ and this->ptemp_rho_.
 */
void ESolver_OF::allocate_array()
{
    // Initialize the "wavefunction", which is sqrt(rho)
    this->psi_ = new psi::Psi<double>(1, GlobalV::NSPIN, this->pw_rho->nrxx);
    ModuleBase::Memory::record("OFDFT::Psi", sizeof(double) * GlobalV::NSPIN * this->pw_rho->nrxx);
    this->pphi_ = new double*[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        this->pphi_[is] = this->psi_->get_pointer(is);
    }
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT PHI");

    // initialize chemical potential, step length, ...
    delete this->ptemp_rho_;
    this->ptemp_rho_ = new Charge();
    this->ptemp_rho_->set_rhopw(this->pw_rho);
    this->ptemp_rho_->allocate(GlobalV::NSPIN);

    this->mu_ = new double[GlobalV::NSPIN];
    this->theta_ = new double[GlobalV::NSPIN];
    this->pdLdphi_ = new double*[GlobalV::NSPIN];
    this->pdEdphi_ = new double*[GlobalV::NSPIN];
    this->pdirect_ = new double*[GlobalV::NSPIN];
    this->precip_dir_ = new std::complex<double> *[GlobalV::NSPIN];

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        this->pdLdphi_[is] = new double[this->pw_rho->nrxx];
        this->pdEdphi_[is] = new double[this->pw_rho->nrxx];
        this->pdirect_[is] = new double[this->pw_rho->nrxx];
        this->precip_dir_[is] = new std::complex<double>[pw_rho->npw];
    }
    ModuleBase::Memory::record("OFDFT::pdLdphi_", sizeof(double) * GlobalV::NSPIN * this->pw_rho->nrxx);
    ModuleBase::Memory::record("OFDFT::pdEdphi_", sizeof(double) * GlobalV::NSPIN * this->pw_rho->nrxx);
    ModuleBase::Memory::record("OFDFT::pdirect_", sizeof(double) * GlobalV::NSPIN * this->pw_rho->nrxx);
    ModuleBase::Memory::record("OFDFT::precip_dir_", sizeof(std::complex<double>) * GlobalV::NSPIN * this->pw_rho->npw);
}

/**
 * @brief Get dL/dphi = dL/drho * drho/dphi = (dE/drho - mu) * 2 * ptemp_phi and store it in rdLdphi
 * 
 * @param [in] ptemp_phi phi
 * @param [out] rdLdphi dL/dphi
 */
void ESolver_OF::cal_potential(double *ptemp_phi, double *rdLdphi)
{
    double **dEdtemp_phi = new double*[GlobalV::NSPIN];
    double **temp_phi = new double*[GlobalV::NSPIN];

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        dEdtemp_phi[is] = new double[this->pw_rho->nrxx];
        if (is == this->tn_spin_flag_)
        {
            temp_phi[is] = ptemp_phi;
        }
        else
        {
            temp_phi[is] = this->pphi_[is];
        }
        for (int ir = 0; ir < this->pw_rho->nrxx; ++ir)
        {
            this->ptemp_rho_->rho[is][ir] = temp_phi[is][ir] * temp_phi[is][ir];
        }
    }

    if(GlobalV::NSPIN==4) GlobalC::ucell.cal_ux();
    this->pelec->pot->update_from_charge(this->ptemp_rho_, &GlobalC::ucell);
    ModuleBase::matrix& vr_eff = this->pelec->pot->get_effective_v();

    this->kinetic_potential(this->ptemp_rho_->rho, temp_phi, vr_eff);
    for (int i = 0; i < this->pw_rho->nrxx; ++i)
    {
        dEdtemp_phi[this->tn_spin_flag_][i] = vr_eff(this->tn_spin_flag_,i);
    }
    double temp_mu = this->cal_mu(ptemp_phi, dEdtemp_phi[this->tn_spin_flag_], this->nelec_[this->tn_spin_flag_]);
    for (int i = 0; i < this->pw_rho->nrxx; ++i)
    {
        rdLdphi[i] = dEdtemp_phi[this->tn_spin_flag_][i] - 2. * temp_mu * ptemp_phi[i];
    }
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        delete[] dEdtemp_phi[is];
    }
    delete[] dEdtemp_phi;
    delete[] temp_phi;
}

/**
 * @brief Calculate dE/dTheta and store it in rdEdtheta.
 * dE/dTheta = <dE / dtemp_phi | dtemp_phi / dTheta>
 *           = <dE / dtemp_phi | - sin(theta) * phi + cos(theta) * direction>
 * 
 * @param [in] ptemp_phi 
 * @param [in] temp_rho 
 * @param [in] ucell 
 * @param [in] ptheta 
 * @param [out] rdEdtheta dE/dTheta
 */
void ESolver_OF::cal_dEdtheta(double **ptemp_phi, Charge* temp_rho, UnitCell& ucell, double *ptheta, double *rdEdtheta)
{
    double *dphi_dtheta = new double[this->pw_rho->nrxx];

    if(GlobalV::NSPIN==4) ucell.cal_ux();
    this->pelec->pot->update_from_charge(temp_rho, &ucell);
    ModuleBase::matrix& vr_eff = this->pelec->pot->get_effective_v();

    this->kinetic_potential(temp_rho->rho, ptemp_phi, vr_eff);
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        for (int ir = 0; ir < this->pw_rho->nrxx; ++ir)
        {
            this->pdEdphi_[is][ir] = vr_eff(is,ir);
            dphi_dtheta[ir] = - this->pphi_[is][ir] * sin(ptheta[is]) + this->pdirect_[is][ir] * cos(ptheta[is]);
        }
        rdEdtheta[is] = this->inner_product(this->pdEdphi_[is], dphi_dtheta, this->pw_rho->nrxx, this->dV_);
        Parallel_Reduce::reduce_all(rdEdtheta[is]);
    }
    delete[] dphi_dtheta;
}

/**
 * @brief Calculate the chemical potential mu.
 * mu = <dE/dphi|phi> / (2 * nelec)
 * 
 * @param pphi 
 * @param pdEdphi 
 * @param nelec 
 * @return mu
 */
double ESolver_OF::cal_mu(double *pphi, double *pdEdphi, double nelec)
{
    double mu = this->inner_product(pphi, pdEdphi, this->pw_rho->nrxx, this->dV_);
    Parallel_Reduce::reduce_all(mu);
    mu = mu / (2.0*nelec);
    return mu;
}

/**
 * @brief Rotate and renormalize the direction |d>, 
 * make it orthogonal to phi (<d|phi> = 0), and <d|d> = nelec
 */
void ESolver_OF::adjust_direction()
{
    // filter the high frequency term in direction if of_full_pw = false
    if (!GlobalV::of_full_pw)
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            pw_rho->real2recip(this->pdirect_[is], this->precip_dir_[is]);
            pw_rho->recip2real(this->precip_dir_[is], this->pdirect_[is]);
        }
    }

    if (GlobalV::NSPIN == 1)
    {
        double temp_theta = 0; // temp_theta = |d'|/|d0 + phi|, theta = min(theta, temp_theta)

        // (1) make direction orthogonal to phi
        // |d'> = |d0> - |phi><phi|d0>/nelec
        double inner_phi_direction = this->inner_product(this->pphi_[0], this->pdirect_[0], this->pw_rho->nrxx, this->dV_);
        Parallel_Reduce::reduce_all(inner_phi_direction);
        for (int i = 0; i < this->pw_rho->nrxx; ++i)
        {
            temp_theta += pow(this->pdirect_[0][i] + this->pphi_[0][i], 2);
            this->pdirect_[0][i] = this->pdirect_[0][i] - this->pphi_[0][i] * inner_phi_direction / this->nelec_[0];
        }
        Parallel_Reduce::reduce_all(temp_theta);
        temp_theta = std::sqrt(temp_theta);

        // (2) renormalize direction
        // |d> = |d'> * \sqrt(nelec) / <d'|d'>
        double norm_direction = this->inner_product(this->pdirect_[0], this->pdirect_[0], this->pw_rho->nrxx, this->dV_);
        Parallel_Reduce::reduce_all(norm_direction);
        norm_direction = std::sqrt(norm_direction);
        for (int i = 0; i < this->pw_rho->nrxx; ++i)
        {
            this->pdirect_[0][i] = std::sqrt(this->nelec_[0]) * this->pdirect_[0][i] / norm_direction;
        }

        temp_theta = norm_direction/temp_theta;
        this->theta_[0] = std::min(this->theta_[0], temp_theta);
    }
    else if (GlobalV::NSPIN == 2) // theta = 0
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            // (1) make direction orthogonal to phi
            // |d'> = |d0> - |phi><phi|d0>/nelec
            double inner_phi_direction = this->inner_product(this->pphi_[is], this->pdirect_[is], this->pw_rho->nrxx, this->dV_);
            Parallel_Reduce::reduce_all(inner_phi_direction);
            for (int i = 0; i < this->pw_rho->nrxx; ++i)
            {
                this->pdirect_[is][i] = this->pdirect_[is][i] - this->pphi_[is][i] * inner_phi_direction / this->nelec_[is];
            }

            // (2) renormalize direction
            // |d> = |d'> * \sqrt(nelec) / <d'|d'>
            double norm_direction = this->inner_product(this->pdirect_[is], this->pdirect_[is], this->pw_rho->nrxx, this->dV_);
            Parallel_Reduce::reduce_all(norm_direction);
            norm_direction = std::sqrt(norm_direction);
            for (int i = 0; i < this->pw_rho->nrxx; ++i)
            {
                this->pdirect_[is][i] = std::sqrt(this->nelec_[is]) * this->pdirect_[is][i] / norm_direction;
            }
            this->theta_[is] = 0.;
        }
    }
}

/**
 * @brief Make sure that dEdtheta<0 at theta = 0,
 * preparing to call the line search
 * 
 * @param dEdtheta 
 * @param ptemp_phi 
 * @param ucell 
 */
void ESolver_OF::check_direction(double *dEdtheta, double **ptemp_phi, UnitCell& ucell)
{
    double *temp_theta = new double[GlobalV::NSPIN];
    ModuleBase::GlobalFunc::ZEROS(temp_theta, GlobalV::NSPIN);

    double max_dEdtheta = 1e5; // threshould of dEdtheta, avoid the unstable optimization
    this->cal_dEdtheta(ptemp_phi, this->ptemp_rho_, ucell, temp_theta, dEdtheta);

    // Assert dEdtheta(theta = 0) < 0, otherwise line search will not work.
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        if (dEdtheta[is] > max_dEdtheta)
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
            this->adjust_direction();
            this->cal_dEdtheta(ptemp_phi, this->ptemp_rho_, ucell, temp_theta, dEdtheta);
            if (dEdtheta[is] > max_dEdtheta)
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
    delete[] temp_theta;
}

/**
 * @brief ONLY used for test.
 * Check the validity of KEDF
 * 
 * @param dEdtheta 
 * @param ptemp_phi 
 * @param ucell 
 */
void ESolver_OF::test_direction(double *dEdtheta, double **ptemp_phi, UnitCell& ucell)
{
    double temp_energy = 0.;
    if (this->iter_ == 0)
    {
        for (int i = -100; i < 100; ++i)
        {
            this->theta_[0] = 0.001 * i;
            for (int ir = 0; ir < this->pw_rho->nrxx; ++ir)
            {
                ptemp_phi[0][ir] = this->pphi_[0][ir] * cos(this->theta_[0]) + this->pdirect_[0][ir] *
                sin(this->theta_[0]); ptemp_rho_->rho[0][ir] = ptemp_phi[0][ir] * ptemp_phi[0][ir];
            }
            this->cal_dEdtheta(ptemp_phi, ptemp_rho_, ucell, this->theta_, dEdtheta);
            this->pelec->cal_energies(2);
            temp_energy = this->pelec->f_en.etot;
            double kinetic_energy = 0.;
            double pseudopot_energy = 0.;
            kinetic_energy = this->kinetic_energy();
            pseudopot_energy = this->inner_product(this->pelec->pot->get_fixed_v(), this->ptemp_rho_->rho[0], this->pw_rho->nrxx, this->dV_);
            Parallel_Reduce::reduce_all(pseudopot_energy);
            temp_energy += kinetic_energy + pseudopot_energy;
            GlobalV::ofs_warning << i << "    " << dEdtheta[0] << "    " << temp_energy << std::endl;
            if (this->theta_[0] == 0) std::cout << "dEdtheta    " << dEdtheta[0]<< std::endl;
        }
        exit(0);
    }
}

/**
 * @brief Print nessecary information to the screen,
 * and write the components of the total energy into running_log.
 */
void ESolver_OF::print_info()
{
    if (this->iter_ == 0){
        std::cout << "============================== Running OFDFT ==============================" <<  std::endl;
        std::cout << "Iter        Etot(Ha)          Mu(Ha)      Theta      PotNorm     deltaE(Ha)" << std::endl;
        // cout << "======================================== Running OFDFT ========================================" <<  endl;
        // cout << "Iter        Etot(Ha)          Theta       PotNorm        min/max(den)          min/max(dE/dPhi)" << endl;
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
    << std::setw(12) << std::setprecision(3) << this->mu_[0]/2.
    << std::setw(12) << this->theta_[0]
    << std::setw(12) << this->normdLdphi_
    << std::setw(12) << (this->energy_current_ - this->energy_last_)/2. << std::endl;
    // ============ used to compare with PROFESS3.0 ================
    // << setw(10) << minDen << "/ " << setw(12) << maxDen
    // << setw(10) << minPot << "/ " << setw(10) << maxPot << endl;
    // =============================================================

    GlobalV::ofs_running << std::setprecision(12);
    GlobalV::ofs_running << std::setiosflags(std::ios::right);

    GlobalV::ofs_running << "\nIter" << this->iter_ << ": the norm of potential is " << this->normdLdphi_ << std::endl;

    std::vector<std::string> titles;
    std::vector<double> energies_Ry;
    std::vector<double> energies_eV;
    context.set_context({"title", "energy", "energy"});
    if (INPUT.printe > 0 && ((this->iter_ + 1) % INPUT.printe == 0 || this->conv_ || this->iter_ == GlobalV::SCF_NMAX))
    {
        titles.push_back("E_Total");        energies_Ry.push_back(this->pelec->f_en.etot);
        titles.push_back("E_Kinetic");      energies_Ry.push_back(this->pelec->f_en.ekinetic);
        titles.push_back("E_Hartree");      energies_Ry.push_back(this->pelec->f_en.hartree_energy);
        titles.push_back("E_xc");           energies_Ry.push_back(this->pelec->f_en.etxc - this->pelec->f_en.etxcc);
        titles.push_back("E_IonElec");      energies_Ry.push_back(this->pelec->f_en.eion_elec);
        titles.push_back("E_Ewald");        energies_Ry.push_back(this->pelec->f_en.ewald_energy);
        if (this->of_kinetic_ == "tf" || this->of_kinetic_ == "tf+" || this->of_kinetic_ == "wt")
        {
            titles.push_back("TF KEDF"); energies_Ry.push_back(this->tf_->tf_energy);
        }
        if (this->of_kinetic_ == "vw" || this->of_kinetic_ == "tf+" || this->of_kinetic_ == "wt"
            || this->of_kinetic_ == "lkt")
        {
            titles.push_back("vW KEDF"); energies_Ry.push_back(this->vw_->vw_energy);
        }
        if (this->of_kinetic_ == "wt")
        {
            titles.push_back("WT KEDF"); energies_Ry.push_back(this->wt_->wt_energy);
        }
        if (this->of_kinetic_ == "lkt")
        {
            titles.push_back("LKT KEDF"); energies_Ry.push_back(this->lkt_->lkt_energy);
        }
        std::string vdw_method = INPUT.vdw_method;
        if (vdw_method == "d2") // Peize Lin add 2014-04, update 2021-03-09
        {
            titles.push_back("E_vdwD2"); energies_Ry.push_back(this->pelec->f_en.evdw);
        }
        else if (vdw_method == "d3_0" || vdw_method == "d3_bj") // jiyy add 2019-05, update 2021-05-02
        {
            titles.push_back("E_vdwD3"); energies_Ry.push_back(this->pelec->f_en.evdw);
        }
        if (GlobalV::imp_sol)
        {
            titles.push_back("E_sol_el"); energies_Ry.push_back(this->pelec->f_en.esol_el);
            titles.push_back("E_sol_cav"); energies_Ry.push_back(this->pelec->f_en.esol_cav);
        }
        if (GlobalV::EFIELD_FLAG)
        {
            titles.push_back("E_efield"); energies_Ry.push_back(elecstate::Efield::etotefield);
        }
        if (GlobalV::GATE_FLAG)
        {
            titles.push_back("E_gatefield"); energies_Ry.push_back(elecstate::Gatefield::etotgatefield);
        }
    }
    else
    {
        titles.push_back("E_Total"); energies_Ry.push_back(this->pelec->f_en.etot);
    }

    if (GlobalV::TWO_EFERMI)
    {
        titles.push_back("E_Fermi_up"); energies_Ry.push_back(this->mu_[0]);
        titles.push_back("E_Fermi_dw"); energies_Ry.push_back(this->mu_[1]);
    }
    else
    {
        titles.push_back("E_Fermi"); energies_Ry.push_back(this->mu_[0]);
    }
    for (int i = 0; i < titles.size(); ++i)
    {
        energies_eV.push_back(energies_Ry[i] * ModuleBase::Ry_to_eV);
    }
    context.enable_title();
    context << "Energy" << titles << "Rydberg" << energies_Ry << "eV" << energies_eV;
    context.center_title();
    GlobalV::ofs_running << context.str() << std::endl;
}
}