#include "kedf_manager.h"

/**
 * @brief [Interface to kedf]
 * Initialize the KEDFs.
 *
 * @param inp
 * @param pw_rho pw basis for charge density
 * @param dV volume of one grid point in real space
 * @param nelec number of electrons in the unit cell
 */
void KEDF_Manager::init(
    const Input_para& inp,
    ModulePW::PW_Basis* pw_rho,
    const double dV,
    const double nelec
)
{
    this->of_kinetic_ = inp.of_kinetic;

    //! Thomas-Fermi (TF) KEDF, TF+ KEDF, and Want-Teter (WT) KEDF
    if (this->of_kinetic_ == "tf"
     || this->of_kinetic_ == "tf+"
     || this->of_kinetic_ == "wt"
     || this->of_kinetic_ == "ml")
    {
        if (this->tf_ == nullptr)
        {
            this->tf_ = new KEDF_TF();
        }
        this->tf_->set_para(pw_rho->nrxx, dV, inp.of_tf_weight);
    }

    //! vW, TF+, WT, and LKT KEDFs
    if (this->of_kinetic_ == "vw" || this->of_kinetic_ == "tf+" || this->of_kinetic_ == "wt"
        || this->of_kinetic_ == "lkt" || this->of_kinetic_ == "ml")
    {
        if (this->vw_ == nullptr)
        {
            this->vw_ = new KEDF_vW();
        }
        this->vw_->set_para(dV, inp.of_vw_weight);
    }

    //! Wang-Teter KEDF
    if (this->of_kinetic_ == "wt")
    {
        if (this->wt_ == nullptr)
        {
            this->wt_ = new KEDF_WT();
        }
        this->wt_->set_para(dV,
                            inp.of_wt_alpha,
                            inp.of_wt_beta,
                            nelec,
                            inp.of_tf_weight,
                            inp.of_vw_weight,
                            inp.of_wt_rho0,
                            inp.of_hold_rho0,
                            inp.of_read_kernel,
                            inp.of_kernel_file,
                            pw_rho);
    }

    //! LKT KEDF
    if (this->of_kinetic_ == "lkt")
    {
        if (this->lkt_ == nullptr)
        {
            this->lkt_ = new KEDF_LKT();
        }
        this->lkt_->set_para(dV, inp.of_lkt_a);
    }
#ifdef __MLALGO
    if (this->of_kinetic_ == "ml")
    {
        if (this->ml_ == nullptr)
            this->ml_ = new KEDF_ML();
        this->ml_->set_para(pw_rho->nrxx, dV, nelec, inp.of_tf_weight, inp.of_vw_weight, 
                        inp.of_ml_chi_p, inp.of_ml_chi_q, inp.of_ml_chi_xi, inp.of_ml_chi_pnl, inp.of_ml_chi_qnl,
                        inp.of_ml_nkernel, inp.of_ml_kernel, inp.of_ml_kernel_scaling,
                        inp.of_ml_yukawa_alpha, inp.of_ml_kernel_file, inp.of_ml_gamma, inp.of_ml_p, inp.of_ml_q, inp.of_ml_tanhp, inp.of_ml_tanhq,
                        inp.of_ml_gammanl, inp.of_ml_pnl, inp.of_ml_qnl, inp.of_ml_xi, inp.of_ml_tanhxi,
                        inp.of_ml_tanhxi_nl, inp.of_ml_tanh_pnl, inp.of_ml_tanh_qnl, inp.of_ml_tanhp_nl, inp.of_ml_tanhq_nl, inp.of_ml_device, pw_rho);
    }
#endif
}

/**
 * @brief [Interface to kedf]
 * Calculated the kinetic potential and plus it to rpot,
 *
 * @param [in] prho charge density
 * @param [in] pphi phi^2 = rho
 * @param [in] pw_rho pw basis for charge density
 * @param [out] rpot rpot => (rpot + kietic potential) * 2 * pphi
 */
void KEDF_Manager::get_potential(
    const double* const* prho,
    const double* const* pphi,
    ModulePW::PW_Basis* pw_rho,
    ModuleBase::matrix& rpot
)
{

#ifdef __MLALGO
    // for ML KEDF test
    if (PARAM.inp.of_ml_local_test) this->ml_->localTest(prho, pw_rho);
#endif

    if (this->of_kinetic_ == "tf" || this->of_kinetic_ == "tf+" || this->of_kinetic_ == "wt")
    {
        this->tf_->tf_potential(prho, rpot);
    }
    if (this->of_kinetic_ == "wt")
    {
        this->wt_->wt_potential(prho, pw_rho, rpot);
    }
    if (this->of_kinetic_ == "lkt")
    {
        this->lkt_->lkt_potential(prho, pw_rho, rpot);
    }
#ifdef __MLALGO
    if (this->of_kinetic_ == "ml")
    {
        this->ml_->ml_potential(prho, pw_rho, rpot);
        this->tf_->get_energy(prho); // temp
    }
#endif

    // Before call vw_potential, change rpot to rpot * 2 * pphi
    for (int is = 0; is < PARAM.inp.nspin; ++is)
    {
        for (int ir = 0; ir < pw_rho->nrxx; ++ir)
        {
            rpot(is, ir) *= 2.0 * pphi[is][ir];
        }
    }

    if (this->of_kinetic_ == "vw" || this->of_kinetic_ == "tf+" || this->of_kinetic_ == "wt"
        || this->of_kinetic_ == "lkt" || this->of_kinetic_ == "ml")
    {
        this->vw_->vw_potential(pphi, pw_rho, rpot);
    }
}

/**
 * @brief [Interface to kedf]
 * Return the kinetic energy
 *
 * @return kinetic energy
 */
double KEDF_Manager::get_energy()
{
    double kinetic_energy = 0.0;

    if (this->of_kinetic_ == "tf" || this->of_kinetic_ == "tf+" || this->of_kinetic_ == "wt")
    {
        kinetic_energy += this->tf_->tf_energy;
    }

    if (this->of_kinetic_ == "vw" || this->of_kinetic_ == "tf+" || this->of_kinetic_ == "wt"
        || this->of_kinetic_ == "lkt" || this->of_kinetic_ == "ml")
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
#ifdef __MLALGO
    if (this->of_kinetic_ == "ml")
    {
        kinetic_energy += this->ml_->ml_energy;
        if (this->ml_->ml_energy >= this->tf_->tf_energy)
        {
            std::cout << "WARNING: ML >= TF" << std::endl;
            std::cout << "ML Term = " << this->ml_->ml_energy << " Ry, TF Term = " << this->tf_->tf_energy << " Ry." << std::endl;
        }
    }
#endif

    return kinetic_energy;
}

/**
 * @brief [Interface to kedf]
 * Calculated the kinetic energy density, ONLY SPIN=1 SUPPORTED
 *
 * @param [in] prho charge density
 * @param [in] pphi phi = sqrt(rho)
 * @param [in] pw_rho pw basis for charge density
 * @param [out] rtau kinetic energy density
 */
void KEDF_Manager::get_energy_density(
    const double* const* prho,
    const double* const* pphi,
    ModulePW::PW_Basis* pw_rho,
    double** rtau
)
{
    for (int ir = 0; ir < pw_rho->nrxx; ++ir)
    {
        rtau[0][ir] = 0.0;
    }

    if (this->of_kinetic_ == "tf" || this->of_kinetic_ == "tf+" || this->of_kinetic_ == "wt")
    {
        this->tf_->tau_tf(prho, rtau[0]);
    }
    if (this->of_kinetic_ == "vw" || this->of_kinetic_ == "tf+" || this->of_kinetic_ == "wt"
        || this->of_kinetic_ == "lkt")
    {
        this->vw_->tau_vw(pphi, pw_rho, rtau[0]);
    }
    if (this->of_kinetic_ == "wt")
    {
        this->wt_->tau_wt(prho, pw_rho, rtau[0]);
    }
    if (this->of_kinetic_ == "lkt")
    {
        this->lkt_->tau_lkt(prho, pw_rho, rtau[0]);
    }
}

/**
 * @brief [Interface to kedf]
 * Calculate the stress of kedf
 * 
 * @param [in] omega Volume of the unit cell
 * @param [in] prho charge density
 * @param [in] pphi phi^2 = rho
 * @param [in] pw_rho pw basis for charge density
 * @param [out] kinetic_stress_
 */
void KEDF_Manager::get_stress(
    const double omega,
    const double* const* prho,
    const double* const* pphi,
    ModulePW::PW_Basis* pw_rho,
    ModuleBase::matrix& kinetic_stress_
)
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
        this->tf_->get_stress(omega);
        kinetic_stress_ += this->tf_->stress;
    }

    if (this->of_kinetic_ == "vw" || this->of_kinetic_ == "tf+" || this->of_kinetic_ == "wt"
        || this->of_kinetic_ == "lkt")
    {
        this->vw_->get_stress(pphi, pw_rho);
        kinetic_stress_ += this->vw_->stress;
    }

    if (this->of_kinetic_ == "wt")
    {
        this->wt_->get_stress(prho, pw_rho, PARAM.inp.of_vw_weight);
        kinetic_stress_ += this->wt_->stress;
    }

    if (this->of_kinetic_ == "lkt")
    {
        this->lkt_->get_stress(prho, pw_rho);
        kinetic_stress_ += this->lkt_->stress;
    }
    if (this->of_kinetic_ == "ml")
    {
        std::cout << "Sorry, the stress of MPN KEDF is not yet supported." << std::endl;
    }
}

void KEDF_Manager::record_energy(
    std::vector<std::string> &titles,
    std::vector<double> &energies_Ry
)
{
    if (this->of_kinetic_ == "tf" || this->of_kinetic_ == "tf+" || this->of_kinetic_ == "wt")
    {
        titles.push_back("TF KEDF");
        energies_Ry.push_back(this->tf_->tf_energy);
    }
    if (this->of_kinetic_ == "vw" || this->of_kinetic_ == "tf+" || this->of_kinetic_ == "wt"
        || this->of_kinetic_ == "lkt" || this->of_kinetic_ == "ml")
    {
        titles.push_back("vW KEDF");
        energies_Ry.push_back(this->vw_->vw_energy);
    }
    if (this->of_kinetic_ == "wt")
    {
        titles.push_back("WT KEDF");
        energies_Ry.push_back(this->wt_->wt_energy);
    }
    if (this->of_kinetic_ == "lkt")
    {
        titles.push_back("LKT KEDF");
        energies_Ry.push_back(this->lkt_->lkt_energy);
    }
#ifdef __MLALGO
    if (this->of_kinetic_ == "ml")
    {
        titles.push_back("MPN KEDF");
        energies_Ry.push_back(this->ml_->ml_energy);
    }
#endif
}

// In future, this function should be extended to other KEDFs.
void KEDF_Manager::generate_ml_target(
    const double * const *prho,
    ModulePW::PW_Basis *pw_rho,
    const double *veff
)
{
#ifdef __MLALGO
    this->ml_->generateTrainData(prho, pw_rho, veff);
#endif
}
