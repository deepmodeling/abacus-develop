#include "esolver_of.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace ModuleESolver
{

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

//
// Get dL/dphi = dL/drho * drho/dphi = (dE/drho - mu) * 2 * ptemp_phi and store it in rdLdphi
//
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

//
// Calculate dE/dTheta
// dE/dTheta = <dE/dtempPhi|dtempPhi/dTheta>
//           = <dE/dtempPhi|-phi*sin(theta)+d*cos(theta)>
//
void ESolver_OF::cal_dEdtheta(double **ptemp_phi, Charge* tempRho, double *ptheta, double *rdEdtheta)
{
    double *dphi_dtheta = new double[this->pw_rho->nrxx];

    if(GlobalV::NSPIN==4) GlobalC::ucell.cal_ux();
    this->pelec->pot->update_from_charge(tempRho, &GlobalC::ucell);
    ModuleBase::matrix& vr_eff = this->pelec->pot->get_effective_v();

    this->kinetic_potential(tempRho->rho, ptemp_phi, vr_eff);
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

//
// Rotate and renormalize the direction |d>, make it orthogonal to phi, and <d|d> = nelec
//
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

void ESolver_OF::check_direction(double *dEdtheta, double **ptemp_phi)
{
    double *temp_theta = new double[GlobalV::NSPIN];
    ModuleBase::GlobalFunc::ZEROS(temp_theta, GlobalV::NSPIN);

    double max_dEdtheta = 1e5; // threshould of dEdtheta, avoid the unstable optimization
    this->cal_dEdtheta(ptemp_phi, this->ptemp_rho_, temp_theta, dEdtheta);

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
            this->cal_dEdtheta(ptemp_phi, this->ptemp_rho_, temp_theta, dEdtheta);
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

void ESolver_OF::test_direction(double *dEdtheta, double **ptemp_phi)
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
            this->cal_dEdtheta(ptemp_phi, ptemp_rho_, this->theta_, dEdtheta);
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
}