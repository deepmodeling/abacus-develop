#ifndef STOCHASTIC_WF_H
#define STOCHASTIC_WF_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/complexmatrix.h"
#include "sto_che.h"
//qianrui 2021-2-4

//----------------------------------------------
// Generate stochastic wave functions 
//----------------------------------------------

class Stochastic_WF 
{
	public:

    Stochastic_WF();

    ~Stochastic_WF();


	void init(const int & nchi_in,
		const int & nche_sto_in,
		const int & seed_sto_in,
		const double & emax_sto_in,
		const double & emin_sto_in,
		const std::string & stotype_in);

	// ModuleBase::ComplexMatrix may not be a best filetype to store the electronic wave functions
    ModuleBase::ComplexMatrix* chi0;  	// origin stochastic wavefunctions in real space
	ModuleBase::ComplexMatrix* chiortho;	// stochastic wavefunctions after in reciprocal space orthogonalized with KS wavefunctions
	Stochastic_Chebychev sto_che;
	int nchi; 				// Total number of stochatic obitals; unit in a_0^(3/2)
	int nchip; 				// The number of stochatic obitals in current process 
	int nche_sto;			// number of orders for Chebyshev expansion
	int seed_sto;  //random seed 
	std::string stotype;

	int nbands_diag; // number of bands obtained from diagonalization

	int nbands_total; // number of bands in total, nbands_total=nchi+nbands_diag;

	double emax_sto;
	double emin_sto;

	protected:


	// check the properties of chi
    void check_chi(const ModuleBase::ComplexMatrix *chi)const;

};
#endif 
