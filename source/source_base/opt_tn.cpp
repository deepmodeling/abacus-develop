#include "opt_tn.h"

namespace ModuleBase
{

Opt_TN::Opt_TN()
{
    this->mach_prec_ = std::numeric_limits<double>::epsilon(); // get machine precise
}

void Opt_TN::allocate(int nx)
{
    this->nx_ = nx;
    this->cg_.allocate(this->nx_);
}

void Opt_TN::set_para(double dV)
{
    this->dV_ = dV;
    this->cg_.set_para(this->dV_);
}

void Opt_TN::refresh(int nx_new)
{
    this->iter_ = 0;
    if (nx_new != 0)
    {
        this->nx_ = nx_new;
    }
    this->cg_.refresh(nx_new);
}

double Opt_TN::inner_product(double* pa, double* pb, int length)
{
    double innerproduct = BlasConnector::dot(length, pa, 1, pb, 1);
    innerproduct *= this->dV_;
    return innerproduct;
}

double Opt_TN::get_epsilon(double* px, double* pcg_direction)
{
    double epsilon = 0.;
    double xx = this->inner_product(px, px, this->nx_);
    Parallel_Reduce::reduce_all(xx);
    double dd = this->inner_product(pcg_direction, pcg_direction, this->nx_);
    Parallel_Reduce::reduce_all(dd);
    epsilon = 2 * sqrt(this->mach_prec_) * (1 + sqrt(xx)) / sqrt(dd);
    // epsilon = 2 * sqrt(this->mach_prec_) * (1 + sqrt(this->inner_product(px, px, this->nx_)))
    //         / sqrt(this->inner_product(pcg_direction, pcg_direction, this->nx_));
    return epsilon;
}

} // namespace ModuleBase
