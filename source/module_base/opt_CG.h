#ifndef OPT_CG_H
#define OPT_CG_H

#include <iostream>
#include <math.h>
#include "global_function.h"
#include "module_base/parallel_reduce.h"

namespace ModuleBase
{
// 
// A class designed to deal with optimization problems with CG method.
// Three forms of CG methods have been implemented, including standard flow to solve 
// the linear equation Ax = b, Polak-Ribire (PR) form and Hager-Zhang (HZ) form to
// solve general optimization problems min{f(x)}.
// We adopt following abbreviation
// x -> solution
// d -> direction
// g -> gradient
// 
class Opt_CG
{
public:
    Opt_CG();
    ~Opt_CG();

    void init_b(
        double *pinp_b // b in the linear equation Ax = b
    );
    void allocate(
        int nx // length of the solution array x
    );
    void set_para(
        double dV
    );
    void refresh(
        int nx_new=0, // length of new x, default 0 means the length doesn't change
        double *pinp_b=nullptr // new b in Ax = b, default nullptr means we are dealing with general case
    );

    void next_direct(
        double *pgradient, // Ad for linear equaiont Ax=b, and gradient for general case 
        int label, // 0 for solve Ax=b, 1 for PR form, 2 for HZ form
        double *rdirect // next direct
    );
    double step_length(
        double *pAd, // Ad for Ax=b
        double *pdirect, // direct
        int &ifPD // if postive definit
    );

    double get_residual() {return sqrt(this->gg_);};
    int get_iter() {return this->iter_;}

    // void ZEROS(double *x, int n)
    // {
    //     for (int i = 0; i < n; ++i) x[i] =0;
    // }

private:
    double dV_ = 1.;
    int nx_ = 0; // length of the sulotion array x
    int iter_ = 0; // number of iteration
    double gg_ = 1000; // gradient dot gradient
    double beta_ = 0.; // d = -g + beta * d
    double eta_ = 0.01; // a constand used in HZ form
    double *pdirect_old_ = nullptr; // direction of last step
    double *pgradient_old_ = nullptr; // gradient, for meth=0, gradient is minus residual r.

    // only for standard CG
    double alpha_ = 0.; // step length in standard CG
    double *pb_ = nullptr; // b in Ax=b, only for standard CG

    void stantard_CGdirect(
        double *pAd, // Ad for Ax=b
        double *rdirect // next direct
    );
    void PR_beta(
        double *pgradient // df(x)/dx
    );
    void HZ_beta(
        double *pgradient // df(x)/dx
    );
    double inner_product(double *pa, double *pb, int length)
    {
        double innerproduct = BlasConnector::dot(length, pa, 1, pb, 1);
        innerproduct *= this->dV_;
        return innerproduct;
    }
};
}
#endif