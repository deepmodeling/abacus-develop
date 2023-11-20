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
void ESolver_OF::kineticPotential(double** prho, double** pphiInpt, ModuleBase::matrix& rpot)
{
    if (this->of_kinetic_ == "tf" || this->of_kinetic_ == "tf+" || this->of_kinetic_ == "wt")
    {
        this->tf_->tf_potential(prho, rpot);
    }
    if (this->of_kinetic_ == "wt")
    {
        this->wt_->WT_potential(prho, this->pw_rho, rpot);
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
        this->vw_->vW_potential(pphiInpt, this->pw_rho, rpot);
    }
}

// Return the this->of_kinetic_ energy
double ESolver_OF::kineticEnergy()
{
    double kinetic_energy = 0.;

    if (this->of_kinetic_ == "tf" || this->of_kinetic_ == "tf+" || this->of_kinetic_ == "wt")
    {
        kinetic_energy += this->tf_->TFenergy;
    }
    if (this->of_kinetic_ == "vw" || this->of_kinetic_ == "tf+" || this->of_kinetic_ == "wt"
        || this->of_kinetic_ == "lkt")
    {
        kinetic_energy += this->vw_->vWenergy;
    }
    if (this->of_kinetic_ == "wt")
    {
        kinetic_energy += this->wt_->WTenergy;
    }
    if (this->of_kinetic_ == "lkt")
    {
        kinetic_energy += this->lkt_->LKTenergy;
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
} // namespace ModuleESolver