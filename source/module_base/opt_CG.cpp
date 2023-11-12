#include "opt_CG.h"

namespace ModuleBase
{
Opt_CG::Opt_CG()
{
    this->pb = NULL;
    this->pdirect_old = NULL;
    this->pgradient_old = NULL;
}

Opt_CG::~Opt_CG()
{
    if (this->pb != NULL) delete[] this->pb;
    if (this->pdirect_old != NULL) delete[] this->pdirect_old;
    if (this->pgradient_old != NULL) delete[] this->pgradient_old;
}

// 
// Initialize b before solving Ax = b. 
// 
void Opt_CG::init_b(
    double *pinp_b // b in the linear equation Ax = b
)
{
    if (this->pb != NULL) delete[] this->pb;
    this->pb = new double[this->nx];
    for (int i = 0; i < nx; ++i) this->pb[i] = pinp_b[i];
}

// 
// Allocate space for pdirect_old and pgradient_old.
// 
void Opt_CG::allocate(
    int nx // length of the solution array x
)
{
    this->nx = nx;
    if (this->pdirect_old != NULL) delete[] this->pdirect_old;
    if (this->pgradient_old != NULL) delete[] this->pgradient_old;
    this->pdirect_old = new double[this->nx];
    this->pgradient_old = new double[this->nx];
    ModuleBase::GlobalFunc::ZEROS(this->pdirect_old, this->nx);
    ModuleBase::GlobalFunc::ZEROS(this->pgradient_old, this->nx);
}

void Opt_CG::setPara(
    double dV
)
{
    this->dV = dV;
}

// 
// Refresh the class. 
// If nx changes, reallocate space. If b is provided, initialize it.
// 
void Opt_CG::refresh(
    int nx_new, // length of new x, default 0 means the length doesn't change
    double *pinp_b // new b in Ax = b, default NULL means we are dealing with general case
)
{
    this->iter = 0;
    this->alpha = 0.;
    this->beta = 0.;
    if (nx_new!=0)
    {
        this->allocate(nx_new);
    }
    else
    {
        ModuleBase::GlobalFunc::ZEROS(this->pdirect_old, this->nx);
        ModuleBase::GlobalFunc::ZEROS(this->pgradient_old, this->nx);
    }
    if (pinp_b != NULL) this->init_b(pinp_b);
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
        if (this->iter == 0) // if iter == 0, d = -g
        {
            for (int i = 0; i < this->nx; ++i)
            {
                rdirect[i] = - pgradient[i];
                this->pgradient_old[i] = pgradient[i];
                this->pdirect_old[i] = rdirect[i];
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
            for (int i = 0; i < this->nx; ++i)
            {
                rdirect[i] = -pgradient[i] + this->beta * this->pdirect_old[i];
                this->pgradient_old[i] = pgradient[i];
                this->pdirect_old[i] = rdirect[i];
            }
        }
        this->iter++;
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
    double dAd = this->inner_product(pdirect, pAd, this->nx);
    Parallel_Reduce::reduce_all(dAd);
    ifPD = 0;
    // check for positive-definiteness, very important for convergence
    if (dAd == 0)
    {
        this->alpha = 0;
        return 0;
    }
    else if (dAd < 0)
    {
        if (this->iter == 1)
        {
            ifPD = -1;
        }
        else
        {
            ifPD = -2;
        }
    }
    this->alpha = this->gg / dAd;
    return this->alpha;
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
    if (this->iter == 0)
    {
        for (int i = 0; i < this->nx; ++i)
        {   
            this->pgradient_old[i] = - this->pb[i];
            rdirect[i] = this->pb[i];
            this->pdirect_old[i] = this->pb[i];
        }
    }
    else
    {
        double *temp_gradient = new double[this->nx];
        for (int i = 0; i < this->nx; ++i)
        {
            temp_gradient[i] = this->pgradient_old[i] + this->alpha * pAd[i];
        }
        this->beta = this->inner_product(temp_gradient, temp_gradient, this->nx) / this->gg;
        Parallel_Reduce::reduce_all(this->beta);
        for (int i = 0; i < this->nx; ++i)
        {
            this->pgradient_old[i] = temp_gradient[i];
            rdirect[i] =  - this->pgradient_old[i] + this->beta * this->pdirect_old[i];
            this->pdirect_old[i] = rdirect[i];
        }
        delete[] temp_gradient;
    }
    this->gg = this->inner_product(this->pgradient_old, this->pgradient_old, this->nx);
    Parallel_Reduce::reduce_all(this->gg);
    this->iter++;
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
    temp_beta = this->inner_product(pgradient, pgradient, this->nx);
    temp_beta -= this->inner_product(pgradient, this->pgradient_old, this->nx);
    Parallel_Reduce::reduce_all(temp_beta);
    double gg_old = this->inner_product(this->pgradient_old, this->pgradient_old, this->nx);
    Parallel_Reduce::reduce_all(gg_old);
    // temp_beta /= this->inner_product(this->pgradient_old, this->pgradient_old, this->nx);
    temp_beta /= gg_old;
    this->beta = std::max(0., temp_beta);
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
    double *y = new double[this->nx];
    for (int i = 0; i < this->nx; ++i) y[i] = pgradient[i] - this->pgradient_old[i];
    
    double py = this->inner_product(this->pdirect_old, y, this->nx);
    Parallel_Reduce::reduce_all(py);
    double yy = this->inner_product(y, y, this->nx);
    Parallel_Reduce::reduce_all(yy);
    double pg = this->inner_product(this->pdirect_old, pgradient, this->nx);
    Parallel_Reduce::reduce_all(pg);
    double yg = this->inner_product(y, pgradient, this->nx);
    Parallel_Reduce::reduce_all(yg);
    double temp_beta = (yg - 2 * pg * yy / py) /py;

    double pp = this->inner_product(this->pdirect_old, this->pdirect_old, this->nx);
    Parallel_Reduce::reduce_all(pp);
    double gg = this->inner_product(this->pgradient_old, this->pgradient_old, this->nx);
    Parallel_Reduce::reduce_all(gg);
    double temp_eta = -1 / (sqrt(pp) * std::min(this->eta, sqrt(gg)));

    this->beta = std::max(temp_beta, temp_eta);

    delete[] y;
}
}