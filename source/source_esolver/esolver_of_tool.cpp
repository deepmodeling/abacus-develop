#include "esolver_of.h"
#include "source_base/formatter.h"
#include "source_base/memory.h"
#include "source_estate/module_pot/efield.h"
#include "source_estate/module_pot/gatefield.h"
#include "source_io/module_parameter/parameter.h"
#include "source_estate/cal_ux.h"

namespace ModuleESolver
{

/**
 * @brief Initialize this->pelec, as well as this->pelec->pot
 *
 * @param ucell
 */
void ESolver_OF::init_elecstate(UnitCell& ucell)
{
    if (this->pelec == nullptr)
    {
        this->pelec = new elecstate::ElecState((Charge*)(&chr), this->pw_rho, pw_big);
    }

    delete this->pelec->pot;
    this->pelec->pot = new elecstate::Potential(this->pw_rhod,
                                                this->pw_rho,
                                                &ucell,
                                                &(this->locpp.vloc),
                                                &(this->sf),
                                                &(this->solvent),
                                                &(this->pelec->f_en.etxc),
                                                &(this->pelec->f_en.vtxc));
    // There is no Operator in ESolver_OF, register Potentials here!
    std::vector<std::string> pot_register_in;
    if (PARAM.inp.vion_in_h)
    {
        pot_register_in.push_back("local");
    }
    if (PARAM.inp.vh_in_h)
    {
        pot_register_in.push_back("hartree");
    }
    // no variable can choose xc, maybe it is necessary
    pot_register_in.push_back("xc");
    if (PARAM.inp.imp_sol)
    {
        pot_register_in.push_back("surchem");
    }
    if (PARAM.inp.efield_flag)
    {
        pot_register_in.push_back("efield");
    }
    if (PARAM.inp.gate_flag)
    {
        pot_register_in.push_back("gatefield");
    }
    if (PARAM.inp.ml_exx)
    {
        pot_register_in.push_back("ml_exx");
    }
    // only Potential is not empty, Veff and Meta are available
    if (pot_register_in.size() > 0)
    {
        // register Potential by gathered operator
        this->pelec->pot->pot_register(pot_register_in);
    }
}

/**
 * @brief Allocate the arrays, as well as this->psi_ and this->ptemp_rho_.
 */
void ESolver_OF::allocate_array()
{
    // Initialize the "wavefunction", which is sqrt(rho)
    this->psi_ = new psi::Psi<double>(1, 
                                      PARAM.inp.nspin, 
                                      this->pw_rho->nrxx,
                                      this->pw_rho->nrxx,
                                      true);
    ModuleBase::Memory::record("OFDFT::Psi", sizeof(double) * PARAM.inp.nspin * this->pw_rho->nrxx);
    this->pphi_ = new double*[PARAM.inp.nspin];
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        this->pphi_[is] = this->psi_->get_pointer(is);
    }
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT PHI");

    // initialize chemical potential, step length, ...
    delete this->ptemp_rho_;
    this->ptemp_rho_ = new Charge();
    this->ptemp_rho_->set_rhopw(this->pw_rho);
    const bool kin_den = this->ptemp_rho_->kin_density(); // mohan add 20251202
    this->ptemp_rho_->allocate(PARAM.inp.nspin, kin_den);

    this->theta_ = new double[PARAM.inp.nspin];
    this->pdLdphi_ = new double*[PARAM.inp.nspin];
    this->pdEdphi_ = new double*[PARAM.inp.nspin];
    this->pdirect_ = new double*[PARAM.inp.nspin];
    this->precip_dir_ = new std::complex<double>*[PARAM.inp.nspin];

    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        this->pdLdphi_[is] = new double[this->pw_rho->nrxx];
        this->pdEdphi_[is] = new double[this->pw_rho->nrxx];
        this->pdirect_[is] = new double[this->pw_rho->nrxx];
        this->precip_dir_[is] = new std::complex<double>[pw_rho->npw];
    }
    ModuleBase::Memory::record("OFDFT::pdLdphi_", sizeof(double) * PARAM.inp.nspin * this->pw_rho->nrxx);
    ModuleBase::Memory::record("OFDFT::pdEdphi_", sizeof(double) * PARAM.inp.nspin * this->pw_rho->nrxx);
    ModuleBase::Memory::record("OFDFT::pdirect_", sizeof(double) * PARAM.inp.nspin * this->pw_rho->nrxx);
    ModuleBase::Memory::record("OFDFT::precip_dir_", sizeof(std::complex<double>) * PARAM.inp.nspin * this->pw_rho->npw);
}

/**
 * @brief Get dL/dphi = dL/drho * drho/dphi = (dE/drho - mu) * 2 * ptemp_phi and
 * store it in rdLdphi
 *
 * @param [in] ptemp_phi phi
 * @param [out] rdLdphi dL/dphi
 */
void ESolver_OF::cal_potential(double* ptemp_phi, double* rdLdphi, UnitCell& ucell)
{
    double** dEdtemp_phi = new double*[PARAM.inp.nspin];
    double** temp_phi = new double*[PARAM.inp.nspin];

    for (int is = 0; is < PARAM.inp.nspin; ++is)
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

    elecstate::cal_ux(ucell);
    this->pelec->pot->update_from_charge(this->ptemp_rho_, &ucell);
    ModuleBase::matrix& vr_eff = this->pelec->pot->get_eff_v();

    this->kedf_manager_->get_potential(this->ptemp_rho_->rho,
                                       temp_phi,
                                       this->pw_rho,
                                       vr_eff); // KEDF potential
    for (int i = 0; i < this->pw_rho->nrxx; ++i)
    {
        dEdtemp_phi[this->tn_spin_flag_][i] = vr_eff(this->tn_spin_flag_, i);
    }
    double temp_mu = this->cal_mu(ptemp_phi, dEdtemp_phi[this->tn_spin_flag_], this->nelec_[this->tn_spin_flag_]);
    for (int i = 0; i < this->pw_rho->nrxx; ++i)
    {
        rdLdphi[i] = dEdtemp_phi[this->tn_spin_flag_][i] - 2. * temp_mu * ptemp_phi[i];
    }
    for (int is = 0; is < PARAM.inp.nspin; ++is)
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
void ESolver_OF::cal_dEdtheta(double** ptemp_phi, Charge* temp_rho, UnitCell& ucell, double* ptheta, double* rdEdtheta)
{
    double* dphi_dtheta = new double[this->pw_rho->nrxx];

    elecstate::cal_ux(ucell);
    this->pelec->pot->update_from_charge(temp_rho, &ucell);
    ModuleBase::matrix& vr_eff = this->pelec->pot->get_eff_v();

    this->kedf_manager_->get_potential(temp_rho->rho,
                                       ptemp_phi,
                                       this->pw_rho,
                                       vr_eff); // KEDF potential
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        for (int ir = 0; ir < this->pw_rho->nrxx; ++ir)
        {
            this->pdEdphi_[is][ir] = vr_eff(is, ir);
            dphi_dtheta[ir] = -this->pphi_[is][ir] * sin(ptheta[is]) + this->pdirect_[is][ir] * cos(ptheta[is]);
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
double ESolver_OF::cal_mu(double* pphi, double* pdEdphi, double nelec)
{
    double mu = this->inner_product(pphi, pdEdphi, this->pw_rho->nrxx, this->dV_);
    Parallel_Reduce::reduce_all(mu);
    mu = mu / (2.0 * nelec);
    return mu;
}

/**
 * @brief Rotate and renormalize the direction |d>,
 * make it orthogonal to phi (<d|phi> = 0), and <d|d> = nelec
 */
void ESolver_OF::adjust_direction(void)
{
    // filter the high frequency term in direction if of_full_pw = false
    if (!PARAM.inp.of_full_pw)
    {
        for (int is = 0; is < PARAM.inp.nspin; ++is)
        {
            pw_rho->real2recip(this->pdirect_[is], this->precip_dir_[is]);
            pw_rho->recip2real(this->precip_dir_[is], this->pdirect_[is]);
        }
    }

    if (PARAM.inp.nspin == 1)
    {
        double temp_theta = 0; // temp_theta = |d'|/|d0 + phi|, theta = min(theta, temp_theta)

        // (1) make direction orthogonal to phi
        // |d'> = |d0> - |phi><phi|d0>/nelec
        double inner_phi_direction
            = this->inner_product(this->pphi_[0], this->pdirect_[0], this->pw_rho->nrxx, this->dV_);
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
        double norm_direction
            = this->inner_product(this->pdirect_[0], this->pdirect_[0], this->pw_rho->nrxx, this->dV_);
        Parallel_Reduce::reduce_all(norm_direction);
        norm_direction = std::sqrt(norm_direction);
        for (int i = 0; i < this->pw_rho->nrxx; ++i)
        {
            this->pdirect_[0][i] = std::sqrt(this->nelec_[0]) * this->pdirect_[0][i] / norm_direction;
        }

        temp_theta = norm_direction / temp_theta;
        this->theta_[0] = std::min(this->theta_[0], temp_theta);
    }
    else if (PARAM.inp.nspin == 2) // theta = 0
    {
        for (int is = 0; is < PARAM.inp.nspin; ++is)
        {
            // (1) make direction orthogonal to phi
            // |d'> = |d0> - |phi><phi|d0>/nelec
            double inner_phi_direction
                = this->inner_product(this->pphi_[is], this->pdirect_[is], this->pw_rho->nrxx, this->dV_);
            Parallel_Reduce::reduce_all(inner_phi_direction);
            for (int i = 0; i < this->pw_rho->nrxx; ++i)
            {
                this->pdirect_[is][i]
                    = this->pdirect_[is][i] - this->pphi_[is][i] * inner_phi_direction / this->nelec_[is];
            }

            // (2) renormalize direction
            // |d> = |d'> * \sqrt(nelec) / <d'|d'>
            double norm_direction
                = this->inner_product(this->pdirect_[is], this->pdirect_[is], this->pw_rho->nrxx, this->dV_);
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
void ESolver_OF::check_direction(double* dEdtheta, double** ptemp_phi, UnitCell& ucell)
{
    assert(PARAM.inp.nspin > 0);
    double* temp_theta = new double[PARAM.inp.nspin];
    ModuleBase::GlobalFunc::ZEROS(temp_theta, PARAM.inp.nspin);

    double max_dEdtheta = 1e5; // threshould of dEdtheta, avoid the unstable optimization
    this->cal_dEdtheta(ptemp_phi, this->ptemp_rho_, ucell, temp_theta, dEdtheta);

    // Assert dEdtheta(theta = 0) < 0, otherwise line search will not work.
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        if (dEdtheta[is] > max_dEdtheta)
        {
            std::cout << "dEdtheta    " << dEdtheta[is] << std::endl;
            ModuleBase::WARNING_QUIT("esolver_of.cpp", "dE/dtheta is too large.");
        }
        else if (dEdtheta[is] > 0)
        {
            GlobalV::ofs_warning << "ESolver_OF: WARNING "
                                 << "dEdphi > 0, replace direct with steepest descent method." << std::endl;
            for (int ir = 0; ir < this->pw_rho->nrxx; ++ir)
            {
                this->pdirect_[is][ir] = -this->pdLdphi_[is][ir];
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
                GlobalV::ofs_warning << "ESolver_OF: WARNING "
                                     << "when use steepest dencent method, "
                                        "dEdphi > 0, so we might get minimum."
                                     << std::endl;
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
void ESolver_OF::test_direction(double* dEdtheta, double** ptemp_phi, UnitCell& ucell)
{
    double temp_energy = 0.;
    if (this->iter_ == 0)
    {
        for (int i = -100; i < 100; ++i)
        {
            this->theta_[0] = 0.001 * i;
            for (int ir = 0; ir < this->pw_rho->nrxx; ++ir)
            {
                ptemp_phi[0][ir]
                    = this->pphi_[0][ir] * cos(this->theta_[0]) + this->pdirect_[0][ir] * sin(this->theta_[0]);
                ptemp_rho_->rho[0][ir] = ptemp_phi[0][ir] * ptemp_phi[0][ir];
            }
            this->cal_dEdtheta(ptemp_phi, ptemp_rho_, ucell, this->theta_, dEdtheta);
            this->pelec->cal_energies(2);
            temp_energy = this->pelec->f_en.etot;
            double kinetic_energy = 0.;
            double pseudopot_energy = 0.;
            kinetic_energy = this->kedf_manager_->get_energy();
            pseudopot_energy = this->inner_product(this->pelec->pot->get_fixed_v(),
                                                   this->ptemp_rho_->rho[0],
                                                   this->pw_rho->nrxx,
                                                   this->dV_);
            Parallel_Reduce::reduce_all(pseudopot_energy);
            temp_energy += kinetic_energy + pseudopot_energy;
            GlobalV::ofs_warning << i << "    " << dEdtheta[0] << "    " << temp_energy << std::endl;
			if (this->theta_[0] == 0) 
			{
				std::cout << "dEdtheta    " << dEdtheta[0] << std::endl;
			}
		}
        exit(0);
    }
}

} // namespace ModuleESolver
