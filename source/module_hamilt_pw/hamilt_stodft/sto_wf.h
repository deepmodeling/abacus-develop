#ifndef STOCHASTIC_WF_H
#define STOCHASTIC_WF_H

#include "module_cell/klist.h"
#include "module_base/complexmatrix.h"
#include "module_psi/psi.h"

//----------------------------------------------
// Generate stochastic wave functions 
//----------------------------------------------
class Stochastic_WF 
{
	public:

    Stochastic_WF();

    ~Stochastic_WF();


	void init(K_Vectors* p_kv, const int npwx_in);

    // origin stochastic wavefunctions in real space
    psi::Psi<std::complex<double>>* chi0 = nullptr;
    // stochastic wavefunctions after in reciprocal space orthogonalized with KS wavefunctions
    psi::Psi<std::complex<double>>* chiortho = nullptr;
    // sqrt(f(H))|chi>
    psi::Psi<std::complex<double>>* shchi = nullptr;
    int nchi = 0; 				// Total number of stochatic obitals
	int *nchip = nullptr; 				// The number of stochatic orbitals in current process of each k point.
	int nchip_max = 0;             // Max number of stochastic orbitals among all k points.
	int nks = 0;           //number of k-points
	int *ngk = nullptr;    // ngk in klist
	int npwx = 0;          // max ngk[ik] in all processors

	int nbands_diag = 0; // number of bands obtained from diagonalization

	int nbands_total = 0; // number of bands in total, nbands_total=nchi+nbands_diag;

};
//init stochastic orbitals
void Init_Sto_Orbitals(Stochastic_WF& stowf, const int seed_in);
//update stochastic orbitals
void Update_Sto_Orbitals(Stochastic_WF& stowf, const int seed_in);
//init complete orbitals
void Init_Com_Orbitals(Stochastic_WF& stowf);
#endif 
