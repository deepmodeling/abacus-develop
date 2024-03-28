#include "esolver_of.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace ModuleESolver
{
/**
 * @brief [Interface to kedf]
 * Initialize the KEDFs.
 *
 * @param inp
 */
void ESolver_OF::init_kedf(Input& inp)
{
    //! Thomas-Fermi (TF) KEDF, TF+ KEDF, and Want-Teter (WT) KEDF
    if (this->of_kinetic_ == "tf" || this->of_kinetic_ == "tf+" || this->of_kinetic_ == "wt")
    {
		if (this->tf_ == nullptr)
		{
			this->tf_ = new KEDF_TF();
		}
        this->tf_->set_para(this->pw_rho->nrxx, this->dV_, inp.of_tf_weight);
    }

    //! vW, TF+, WT, and LKT KEDFs
    if (this->of_kinetic_ == "vw" 
     || this->of_kinetic_ == "tf+" 
     || this->of_kinetic_ == "wt"
     || this->of_kinetic_ == "lkt")
    {
		if (this->vw_ == nullptr)
		{
			this->vw_ = new KEDF_vW();
		}
        this->vw_->set_para(this->dV_, inp.of_vw_weight);
    }

    //! Wang-Teter KEDF
    if (this->of_kinetic_ == "wt")
    {
		if (this->wt_ == nullptr)
		{
			this->wt_ = new KEDF_WT();
		}
        this->wt_->set_para(this->dV_,
                            inp.of_wt_alpha,
                            inp.of_wt_beta,
                            this->nelec_[0],
                            inp.of_tf_weight,
                            inp.of_vw_weight,
                            inp.of_wt_rho0,
                            inp.of_hold_rho0,
                            inp.of_read_kernel,
                            inp.of_kernel_file,
                            this->pw_rho);
    }

    //! LKT KEDF
    if (this->of_kinetic_ == "lkt")
    {
		if (this->lkt_ == nullptr)
		{
			this->lkt_ = new KEDF_LKT();
		}
        this->lkt_->set_para(this->dV_, inp.of_lkt_a);
    }
}

/**
 * @brief [Interface to kedf]
 * Calculated the kinetic potential and plus it to rpot,
 *
 * @param [in] prho charge density
 * @param [in] pphi phi^2 = rho
 * @param [out] rpot rpot => (rpot + kietic potential) * 2 * pphi
 */
void ESolver_OF::kinetic_potential(double** prho, double** pphi, ModuleBase::matrix& rpot)
{
    if (this->of_kinetic_ == "tf" || this->of_kinetic_ == "tf+" || this->of_kinetic_ == "wt")
    {
        this->tf_->tf_potential(prho, rpot);
    }
    if (this->of_kinetic_ == "wt")
    {
        this->wt_->wt_potential(prho, this->pw_rho, rpot);
    }
    if (this->of_kinetic_ == "lkt")
    {
        this->lkt_->lkt_potential(prho, this->pw_rho, rpot);
    }

    // Before call vw_potential, change rpot to rpot * 2 * pphi
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        for (int ir = 0; ir < this->pw_rho->nrxx; ++ir)
        {
            rpot(is, ir) *= 2.0 * pphi[is][ir];
        }
    }

    if (this->of_kinetic_ == "vw" 
     || this->of_kinetic_ == "tf+" 
     || this->of_kinetic_ == "wt"
     || this->of_kinetic_ == "lkt")
    {
        this->vw_->vw_potential(pphi, this->pw_rho, rpot);
    }
}

/**
 * @brief [Interface to kedf]
 * Return the kinetic energy
 *
 * @return kinetic energy
 */
double ESolver_OF::kinetic_energy()
{
    double kinetic_energy = 0.0;

    if (this->of_kinetic_ == "tf" 
     || this->of_kinetic_ == "tf+" 
     || this->of_kinetic_ == "wt")
    {
        kinetic_energy += this->tf_->tf_energy;
    }

    if (this->of_kinetic_ == "vw" 
     || this->of_kinetic_ == "tf+" 
     || this->of_kinetic_ == "wt"
     || this->of_kinetic_ == "lkt")
    {
        kinetic_energy += this->vw_->vw_energy;
    }

    if (this->of_kinetic_ == "wt")
    {
        kinetic_energy += this->wt_->wt_energy;
    }

    if (this->of_kinetic_ == "lkt")
    {
        kinetic_energy += this->lkt_->lkt_energy;
    }

    return kinetic_energy;
}

/**
 * @brief [Interface to kedf]
 * Calculate the stress of kedf
 *
 * @param [out] kinetic_stress_
 */
void ESolver_OF::kinetic_stress(ModuleBase::matrix& kinetic_stress_)
{
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            kinetic_stress_(i, j) = 0.0;
        }
    }

    if (this->of_kinetic_ == "tf" || this->of_kinetic_ == "tf+" || this->of_kinetic_ == "wt")
    {
        this->tf_->get_stress(this->pelec->omega);
        kinetic_stress_ += this->tf_->stress;
    }

    if (this->of_kinetic_ == "vw" || this->of_kinetic_ == "tf+" || this->of_kinetic_ == "wt"
        || this->of_kinetic_ == "lkt")
    {
        this->vw_->get_stress(this->pphi_, this->pw_rho);
        kinetic_stress_ += this->vw_->stress;
    }

    if (this->of_kinetic_ == "wt")
    {
        this->wt_->get_stress(pelec->charge->rho, this->pw_rho, GlobalV::of_vw_weight);
        kinetic_stress_ += this->wt_->stress;
    }

    if (this->of_kinetic_ == "lkt")
    {
        this->lkt_->get_stress(pelec->charge->rho, this->pw_rho);
        kinetic_stress_ += this->lkt_->stress;
    }
}

/**
 * @brief [Interface to opt]
 * Initialize the opts
 */
void ESolver_OF::init_opt()
{
	if (this->opt_dcsrch_ == nullptr)
	{
		this->opt_dcsrch_ = new ModuleBase::Opt_DCsrch();
	}

    if (this->of_method_ == "tn")
    {
		if (this->opt_tn_ == nullptr)
		{
			this->opt_tn_ = new ModuleBase::Opt_TN();
		}
        this->opt_tn_->allocate(this->pw_rho->nrxx);
        this->opt_tn_->set_para(this->dV_);
    }
    else if (this->of_method_ == "cg1" || this->of_method_ == "cg2")
    {
		if (this->opt_cg_ == nullptr)
		{
			this->opt_cg_ = new ModuleBase::Opt_CG();
		}
		this->opt_cg_->allocate(this->pw_rho->nrxx);
        this->opt_cg_->set_para(this->dV_);
        this->opt_dcsrch_->set_paras(1e-4, 1e-2);
    }
    else if (this->of_method_ == "bfgs")
    {
        ModuleBase::WARNING_QUIT("esolver_of", "BFGS is not supported now.");
        return;
    }

    // optimize theta if nspin=2
    if (GlobalV::NSPIN == 2)
    {
        this->opt_cg_mag_ = new ModuleBase::Opt_CG;
        this->opt_cg_mag_->allocate(GlobalV::NSPIN);
    }
}

/**
 * @brief [Interface to opt]
 * Call optimization methods to get the optimization direction
 */
void ESolver_OF::get_direction()
{
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        if (this->of_method_ == "tn")
        {
            this->tn_spin_flag_ = is;
            opt_tn_->next_direct(this->pphi_[is],
                                 this->pdLdphi_[is],
                                 this->flag_,
                                 this->pdirect_[is],
                                 this,
                                 &ESolver_OF::cal_potential);
        }
        else if (this->of_method_ == "cg1")
        {
            opt_cg_->next_direct(this->pdLdphi_[is], 1, this->pdirect_[is]);
        }
        else if (this->of_method_ == "cg2")
        {
            opt_cg_->next_direct(this->pdLdphi_[is], 2, this->pdirect_[is]);
        }
        else if (this->of_method_ == "bfgs")
        {
            return;
        }
        else
        {
            ModuleBase::WARNING_QUIT("ESolver_OF", "of_method must be one of CG, TN, or BFGS.");
        }
    }
}

/**
 * @brief [Interface to opt]
 * Call line search to find the best step length
 *
 * @param dEdtheta d E / d theta
 * @param ptemp_phi
 * @param ucell
 */
void ESolver_OF::get_step_length(double* dEdtheta, double** ptemp_phi, UnitCell& ucell)
{
    double temp_energy = 0.0;      // energy of temp_phi and temp_rho
    double kinetic_energy = 0.0;   // kinetic energy
    double pseudopot_energy = 0.0; // electron-ion interaction energy

    if (GlobalV::NSPIN == 1)
    {
        int numDC = 0; // iteration number of line search
        strcpy(this->task_, "START");
        while (true)
        {
            // update energy
            this->pelec->cal_energies(2);
            temp_energy = this->pelec->f_en.etot;
            kinetic_energy = this->kinetic_energy();
            pseudopot_energy = this->inner_product(this->pelec->pot->get_fixed_v(),
                                                   this->ptemp_rho_->rho[0],
                                                   this->pw_rho->nrxx,
                                                   this->dV_);
            Parallel_Reduce::reduce_all(pseudopot_energy);
            temp_energy += kinetic_energy + pseudopot_energy;

            // line search to update theta[0]
            this->opt_dcsrch_->dcSrch(temp_energy, dEdtheta[0], this->theta_[0], this->task_);
            numDC++;

            // decide what to do next according to the output of line search
            if (strncmp(this->task_, "FG", 2) == 0) // continue line search
            {
                // update tempPhi and tempRho
                for (int i = 0; i < this->pw_rho->nrxx; ++i)
                {
                    ptemp_phi[0][i]
                        = this->pphi_[0][i] * cos(this->theta_[0]) + this->pdirect_[0][i] * sin(this->theta_[0]);
                    this->ptemp_rho_->rho[0][i] = ptemp_phi[0][i] * ptemp_phi[0][i];
                }

                // get dEdtheta of new tempPhi and tempRho
                this->cal_dEdtheta(ptemp_phi, this->ptemp_rho_, ucell, this->theta_, dEdtheta);

                if (numDC > this->max_dcsrch_)
                {
                    GlobalV::ofs_warning << "ESolver_OF linesearch: WARNING "
                                         << "excedd the max iter number." << std::endl;
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
        //     double *temp_theta = new double[GlobalV::NSPIN];
        //     ModuleBase::GlobalFunc::ZEROS(pthetaDir, GlobalV::NSPIN);
        //     ModuleBase::GlobalFunc::ZEROS(temp_theta, GlobalV::NSPIN);
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
        //             temp_energy = this->pelec->f_en.etot;
        //             kinetic_energy = this->kinetic_energy();
        //             pseudopot_energy = 0.;
        //             for (int is = 0; is < GlobalV::NSPIN; ++is) {
        //                 pseudopot_energy += this->inner_product(GlobalC::pot.vltot, ptemp_rho_[is],
        //                 this->pw_rho->nrxx, this->dV_);
        //             }
        //             Parallel_Reduce::reduce_all(pseudopot_energy);
        //             temp_energy += kinetic_energy + pseudopot_energy;
        //             this->opt_dcsrch_->dcSrch(temp_energy, dEdalpha, thetaAlpha, this->task_);
        //             numDC++;

        //             if (strncmp(this->task_, "FG", 2) == 0)
        //             {
        //                 for (int is = 0; is < GlobalV::NSPIN; ++is)
        //                 {
        //                     temp_theta[is] = this->theta_[is] + thetaAlpha * pthetaDir[is];
        //                     for (int ir = 0; ir < this->pw_rho->nrxx; ++ir)
        //                     {
        //                         ptemp_phi[is][ir] = this->pphi_[is][ir] * cos(temp_theta[is]) +
        //                         this->pdirect_[is][ir] * sin(temp_theta[is]); ptemp_rho_[is][ir] = ptemp_phi[is][ir]
        //                         * ptemp_phi[is][ir];
        //                     }
        //                 }
        //                 this->cal_dEdtheta(ptemp_phi, ptemp_rho_, temp_theta, dEdtheta);
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
        //     delete[] temp_theta;
        //     delete[] pthetaDir;
        // ========================== Under testing ==========================
    }
    else if (GlobalV::NSPIN == 4)
    {
        ModuleBase::WARNING_QUIT("esolver_of", "Sorry, SPIN4 case is not supported by OFDFT for now.");
    }
}
} // namespace ModuleESolver
