#ifndef STOCHASTIC_WF_H
#define STOCHASTIC_WF_H

#include "klist.h"
#include "../module_base/complexmatrix.h"

//----------------------------------------------
// Generate stochastic wave functions 
//----------------------------------------------
class Stochastic_WF 
{
	public:

    Stochastic_WF();

    ~Stochastic_WF();


	void init(const int nks);

	// ModuleBase::ComplexMatrix may not be a best filetype to store the electronic wave functions
    ModuleBase::ComplexMatrix* chi0;  	// origin stochastic wavefunctions in real space
	ModuleBase::ComplexMatrix* chiortho;	// stochastic wavefunctions after in reciprocal space orthogonalized with KS wavefunctions
	ModuleBase::ComplexMatrix* shchi;     // sqrt(f(H))|chi>
	int nchi; 				// Total number of stochatic obitals
	int *nchip; 				// The number of stochatic orbitals in current process of each k point.
	int nchip_max = 0;             // Max number of stochastic orbitals among all k points.
	int nks;           //number of k-points

	int nbands_diag; // number of bands obtained from diagonalization

	int nbands_total; // number of bands in total, nbands_total=nchi+nbands_diag;

};
//init stochastic orbitals
void Init_Sto_Orbitals(Stochastic_WF& stowf, const int seed_in);
//update stochastic orbitals
void Update_Sto_Orbitals(Stochastic_WF& stowf, const int seed_in);
//init complete orbitals
void Init_Com_Orbitals(Stochastic_WF& stowf, K_Vectors& kv);
#endif 
