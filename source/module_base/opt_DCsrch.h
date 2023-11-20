#ifndef OPT_DCSRCH_H
#define OPT_DCSRCH_H

#include <iostream>
#include "constants.h"

namespace ModuleBase
{

// 
// A interface of line search
// 
class Opt_DCsrch
{
public:
	Opt_DCsrch()
	{
		this->isave_ = new int[3];
    	this->dsave_ = new double[14];
	}
	~Opt_DCsrch()
	{
		delete[] this->isave_; 
		delete[] this->dsave_; 
	}

	// 
	// Reset following parameters.
	// ftol: nonnegative tolerance for the sufficient decrease condition.
	// gtol: nonnegative tolerance for the curvature condition.
	// xtol: nonnegative relative tolerance for an acceptable step. 
	//       The subroutine exits with a warning if the relative difference between sty and stx is less than xtol.
	// stpmin: nonnegative lower bound for the step.
	// stpmax: nonnegative upper bound for the step.
	// 
	// The default setting in PROFESS is
	// TN   ftol = 1e-4   gtol = 2e-1
    // CG   ftol = 1e-4   gtol = 1e-2
    // BFGS ftol = 1e-4   gtol = 2e-1
	// 
	void set_paras(
		double ftol=1e-4, 
		double gtol=2e-1, 
		double xtol=1e-12, 
		double stpmin=0., 
		double stpmax=ModuleBase::PI
	)
	{
		this->ftol_ = ftol;
		this->gtol_ = gtol;
		this->xtol_ = xtol;
		this->stpmin_ = stpmin;
		this->stpmax_ = stpmax;
	}

	// 
	//  Interface of dcsrch, finding the optimal step length with line search.
	//  Input:
	//  f is the value of the function at x on initial entry.
	//       On subsequent entries f is the value of the function at x + stp * d.
	//  g is the derivative of the function at 0 on initial entry.
	//       On subsequent entries g is the derivative of the function at x + stp * d.
	//  Output:
	//  rstp is the optimized step length, assert the initial value is larger than zero.
	//  rtask is a character variable of length at least 60.
	//        On initial entry task must be set to 'START'. 
	//        On exit task indicates the required action:
	//        If task(1:2) = 'FG' then evaluate the function and derivative at stp and call dcsrch again.
	//  		 If task(1:4) = 'CONV' then the search is successful.
	//        If task(1:4) = 'WARN' then the subroutine is not able to satisfy the convergence conditions.
	//           The exit value of stp contains the best point found during the search.
	//        If task(1:5) = 'ERROR' then there is an error in the input arguments.
	//  
	void dcSrch(
		double &f, 
		double &g, 
		double &rstp, 
		char *rtask
	);

private:
	double ftol_ = 1e-4; // nonnegative tolerance for the sufficient decrease condition.
	double gtol_ = 2e-1; // nonnegative tolerance for the curvature condition.
	double xtol_ = 1e-12; // nonnegative relative tolerance for an acceptable step. The subroutine exits with a warning if the relative difference between sty and stx is less than xtol.
	double stpmin_ = 0.; // nonnegative lower bound for the step.
	double stpmax_ = ModuleBase::PI; // nonnegative upper bound for the step.
	int *isave_ = nullptr; // an integer work array of dimension 2.
    double *dsave_ = nullptr; // a double precision work array of dimension 13.
};
}

#endif