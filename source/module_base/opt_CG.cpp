#include "opt_CG.h"

namespace ModuleBase
{
Opt_CG::Opt_CG(){}

Opt_CG::~Opt_CG()
{
    delete[] this->pb_;
    delete[] this->pdirect_old_;
    delete[] this->pgradient_old_;
}

// 
// Initialize b before solving Ax = b. 
// 
void Opt_CG::init_b(
    double *pinp_b // b in the linear equation Ax = b
)
{
    if (this->pb_ != nullptr) delete[] this->pb_;
    this->pb_ = new double[this->nx_];
    for (int i = 0; i < this->nx_; ++i) this->pb_[i] = pinp_b[i];
}

// 
// Allocate space for pdirect_old and pgradient_old.
// 
void Opt_CG::allocate(
    int nx // length of the solution array x
)
{
    this->nx_ = nx;
    delete[] this->pdirect_old_;
    delete[] this->pgradient_old_;
    this->pdirect_old_ = new double[this->nx_];
    this->pgradient_old_ = new double[this->nx_];
    ModuleBase::GlobalFunc::ZEROS(this->pdirect_old_, this->nx_);
    ModuleBase::GlobalFunc::ZEROS(this->pgradient_old_, this->nx_);
}

void Opt_CG::set_para(
    double dV
)
{
    this->dV_ = dV;
}

// 
// Refresh the class. 
// If nx changes, reallocate space. If b is provided, initialize it.
// 
void Opt_CG::refresh(
    int nx_new, // length of new x, default 0 means the length doesn't change
    double *pinp_b // new b in Ax = b, default nullptr means we are dealing with general case
)
{
    this->iter_ = 0;
    this->alpha_ = 0.;
    this->beta_ = 0.;
    if (nx_new!=0)
    {
        this->allocate(nx_new);
    }
    else
    {
        ModuleBase::GlobalFunc::ZEROS(this->pdirect_old_, this->nx_);
        ModuleBase::GlobalFunc::ZEROS(this->pgradient_old_, this->nx_);
    }
    if (pinp_b != nullptr) this->init_b(pinp_b);
}

// 
// Get next optimization direction.
// Input:
// pgradient: Ad for linear equaiont Ax=b, and gradient for general case 
// label: 0 for solve Ax=b, 1 for PR form, 2 for HZ form.
// 
void Opt_CG::next_direct(
    double *pgradient, // Ad for linear equaiont Ax=b, and gradient for general case 
    int label, // 0 for solve Ax=b, 1 for PR form, 2 for HZ form
    double *rdirect // next direct
)
{
    if (label == 0) // standard CG to solve Ap=x
    {
        this->stantard_CGdirect(pgradient, rdirect);
    }
    else if (label == 1 or label == 2) // FR formula or HZ form
    {
        if (this->iter_ == 0) // if iter == 0, d = -g
        {
            for (int i = 0; i < this->nx_; ++i)
            {
                rdirect[i] = - pgradient[i];
                this->pgradient_old_[i] = pgradient[i];
                this->pdirect_old_[i] = rdirect[i];
            }
        }
        else // d = -g + beta * d
        {
            if (label == 1)
            {
                this->PR_beta(pgradient);
            }
            else if (label == 2)
            {
                this->HZ_beta(pgradient);
            }
            for (int i = 0; i < this->nx_; ++i)
            {
                rdirect[i] = -pgradient[i] + this->beta_ * this->pdirect_old_[i];
                this->pgradient_old_[i] = pgradient[i];
                this->pdirect_old_[i] = rdirect[i];
            }
        }
        this->iter_++;
    }
}

// 
// Get step length, only work for standard CG.
// alpha = rr/dAd
// 
double Opt_CG::step_length(
    double *pAd, // Ad for Ax=b
    double *pdirect, // direct
    int &ifPD // 0 if positive definite, -1, -2 when not
)
{
    double dAd = this->inner_product(pdirect, pAd, this->nx_);
    Parallel_Reduce::reduce_all(dAd);
    ifPD = 0;
    // check for positive-definiteness, very important for convergence
    if (dAd == 0)
    {
        this->alpha_ = 0;
        return 0;
    }
    else if (dAd < 0)
    {
        if (this->iter_ == 1)
        {
            ifPD = -1;
        }
        else
        {
            ifPD = -2;
        }
    }
    this->alpha_ = this->gg_ / dAd;
    return this->alpha_;
}

//
// Get next optimization direction with standard CG workflow.
// Only work for solving Ax=b. 
//
void Opt_CG::stantard_CGdirect(        
    double *pAd, // Ad for Ax=b
    double *rdirect // next direct
)
{
    if (this->iter_ == 0)
    {
        for (int i = 0; i < this->nx_; ++i)
        {   
            this->pgradient_old_[i] = - this->pb_[i];
            rdirect[i] = this->pb_[i];
            this->pdirect_old_[i] = this->pb_[i];
        }
    }
    else
    {
        double *temp_gradient = new double[this->nx_];
        for (int i = 0; i < this->nx_; ++i)
        {
            temp_gradient[i] = this->pgradient_old_[i] + this->alpha_ * pAd[i];
        }
        this->beta_ = this->inner_product(temp_gradient, temp_gradient, this->nx_) / this->gg_;
        Parallel_Reduce::reduce_all(this->beta_);
        for (int i = 0; i < this->nx_; ++i)
        {
            this->pgradient_old_[i] = temp_gradient[i];
            rdirect[i] =  - this->pgradient_old_[i] + this->beta_ * this->pdirect_old_[i];
            this->pdirect_old_[i] = rdirect[i];
        }
        delete[] temp_gradient;
    }
    this->gg_ = this->inner_product(this->pgradient_old_, this->pgradient_old_, this->nx_);
    Parallel_Reduce::reduce_all(this->gg_);
    this->iter_++;
}

// 
// Get beta in PR form.
// beta_k = max{0, <g_k, g_k-g_{k-1}>/<g_{k-1}, g_{k-1}>}
// <> means inner product.
// 
void Opt_CG::PR_beta(
    double *pgradient // df(x)/dx
)
{
    double temp_beta = 0.;
    temp_beta = this->inner_product(pgradient, pgradient, this->nx_);
    temp_beta -= this->inner_product(pgradient, this->pgradient_old_, this->nx_);
    Parallel_Reduce::reduce_all(temp_beta);
    double gg_old = this->inner_product(this->pgradient_old_, this->pgradient_old_, this->nx_);
    Parallel_Reduce::reduce_all(gg_old);
    // temp_beta /= this->inner_product(this->pgradient_old_, this->pgradient_old_, this->nx_);
    temp_beta /= gg_old;
    this->beta_ = std::max(0., temp_beta);
}

// 
// Get beta in HZ form.
// See formula in 
// Hager W W, Zhang H. SIAM Journal on optimization, 2005, 16(1): 170-192
// 
void Opt_CG::HZ_beta(
    double *pgradient // df(x)/dx
)
{
    double *y = new double[this->nx_];
    for (int i = 0; i < this->nx_; ++i) y[i] = pgradient[i] - this->pgradient_old_[i];
    
    double py = this->inner_product(this->pdirect_old_, y, this->nx_);
    Parallel_Reduce::reduce_all(py);
    double yy = this->inner_product(y, y, this->nx_);
    Parallel_Reduce::reduce_all(yy);
    double pg = this->inner_product(this->pdirect_old_, pgradient, this->nx_);
    Parallel_Reduce::reduce_all(pg);
    double yg = this->inner_product(y, pgradient, this->nx_);
    Parallel_Reduce::reduce_all(yg);
    double temp_beta = (yg - 2 * pg * yy / py) /py;

    double pp = this->inner_product(this->pdirect_old_, this->pdirect_old_, this->nx_);
    Parallel_Reduce::reduce_all(pp);
    double gg = this->inner_product(this->pgradient_old_, this->pgradient_old_, this->nx_);
    Parallel_Reduce::reduce_all(gg);
    double temp_eta = -1 / (sqrt(pp) * std::min(this->eta_, sqrt(gg)));

    this->beta_ = std::max(temp_beta, temp_eta);

    delete[] y;
}
}