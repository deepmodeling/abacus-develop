#include "esolver_of.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace ModuleESolver
{
void ESolver_OF::init_kedf()
{
    if (this->of_kinetic_ == "tf" || this->of_kinetic_ == "tf+" || this->of_kinetic_ == "wt")
    {
        if (this->tf_ == nullptr)
            this->tf_ = new KEDF_TF();
        this->tf_->set_para(this->pw_rho->nrxx, this->dV_, GlobalV::of_tf_weight);
    }
    if (this->of_kinetic_ == "vw" || this->of_kinetic_ == "tf+" || this->of_kinetic_ == "wt"
        || this->of_kinetic_ == "lkt")
    {
        if (this->vw_ == nullptr)
            this->vw_ = new KEDF_vW();
        this->vw_->set_para(this->dV_, GlobalV::of_vw_weight);
    }
    if (this->of_kinetic_ == "wt")
    {
        if (this->wt_ == nullptr)
            this->wt_ = new KEDF_WT();
        this->wt_->set_para(this->dV_,
                            GlobalV::of_wt_alpha,
                            GlobalV::of_wt_beta,
                            this->nelec_[0],
                            GlobalV::of_tf_weight,
                            GlobalV::of_vw_weight,
                            GlobalV::of_read_kernel,
                            GlobalV::of_kernel_file,
                            this->pw_rho);
    }
    if (this->of_kinetic_ == "lkt")
    {
        if (this->lkt_ == nullptr)
            this->lkt_ = new KEDF_LKT();
        this->lkt_->set_para(this->dV_, GlobalV::of_lkt_a);
    }

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT KEDF");
}

// Calculated this->of_kinetic_ potential and plus it to &rpot, return (rpot + kietic potential) * 2 * pphiInpt
void ESolver_OF::kinetic_potential(double** prho, double** pphiInpt, ModuleBase::matrix& rpot)
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

    // Before call vw_potential, change rpot to rpot * 2 * pphiInpt
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        for (int ir = 0; ir < this->pw_rho->nrxx; ++ir)
        {
            rpot(is, ir) *= 2.0 * pphiInpt[is][ir];
        }
    }

    if (this->of_kinetic_ == "vw" || this->of_kinetic_ == "tf+" || this->of_kinetic_ == "wt"
        || this->of_kinetic_ == "lkt")
    {
        this->vw_->vw_potential(pphiInpt, this->pw_rho, rpot);
    }
}

// Return the this->of_kinetic_ energy
double ESolver_OF::kinetic_energy()
{
    double kinetic_energy = 0.;

    if (this->of_kinetic_ == "tf" || this->of_kinetic_ == "tf+" || this->of_kinetic_ == "wt")
    {
        kinetic_energy += this->tf_->tf_energy;
    }
    if (this->of_kinetic_ == "vw" || this->of_kinetic_ == "tf+" || this->of_kinetic_ == "wt"
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
        this->tf_->get_stress(GlobalC::ucell.omega);
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
        this->wt_->get_stress(GlobalC::ucell.omega, pelec->charge->rho, this->pw_rho, GlobalV::of_vw_weight);
        kinetic_stress_ += this->wt_->stress;
    }
    if (this->of_kinetic_ == "lkt")
    {
        this->lkt_->get_stress(GlobalC::ucell.omega, pelec->charge->rho, this->pw_rho);
        kinetic_stress_ += this->lkt_->stress;
    }
}

void ESolver_OF::init_opt()
{
    if (this->opt_dcsrch_ == nullptr) this->opt_dcsrch_ = new ModuleBase::Opt_DCsrch();

    if (this->of_method_ == "tn")
    {
        if (this->opt_tn_ == nullptr) this->opt_tn_ = new ModuleBase::Opt_TN();
        this->opt_tn_->allocate(this->pw_rho->nrxx);
        this->opt_tn_->set_para(this->dV_);
    }
    else if (this->of_method_ == "cg1" || this->of_method_ == "cg2")
    {
        if (this->opt_cg_ == nullptr) this->opt_cg_ = new ModuleBase::Opt_CG();
        this->opt_cg_->allocate(this->pw_rho->nrxx);
        this->opt_cg_->set_para(this->dV_);
        this->opt_dcsrch_->set_paras(1e-4,1e-2);
    }
    else if (this->of_method_ == "bfgs")
    {
        ModuleBase::WARNING_QUIT("esolver_of", "BFGS is not supported now.");
        return;
    }
}

void ESolver_OF::get_direction()
{
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        if (this->of_method_ == "tn")
        {
            this->tn_spin_flag_ = is;
            opt_tn_->next_direct(this->pphi_[is], this->pdLdphi_[is], this->flag_, this->pdirect_[is], this, &ESolver_OF::calV);
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
} // namespace ModuleESolver