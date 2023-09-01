#ifndef STO_ITER_H
#define STO_ITER_H
#include "sto_wf.h"
#include "module_base/math_chebyshev.h"
#include "sto_hchi.h"
#include "sto_func.h"
#include "module_psi/psi.h"
#include "module_elecstate/elecstate.h"
#include "module_hamilt_general/hamilt.h"

//----------------------------------------------
// Solve for the new electron density and iterate 
// until SCF loop converges (mu, epsilon, rho)
// mu: chemical potential
// epsilon: eigenvalues
// rho: charge density
//----------------------------------------------

class Stochastic_Iter
{

	public:

    // constructor and deconstructor
    Stochastic_Iter();
    ~Stochastic_Iter();

    /**
     * @brief init for iteration process of SDFT
     * 
     * @param method_in 1: slow   2: fast but cost much memories
     * @param pkv_in K_Vectors
     * @param wfc_basis wfc pw basis
     * @param stowf stochastic wave function
     */
    void init(const int method_in, K_Vectors* pkv_in, ModulePW::PW_Basis_K* wfc_basis, Stochastic_WF& stowf);

    void sum_stoband(Stochastic_WF& stowf,
                     elecstate::ElecState* pes,
                     hamilt::Hamilt<double>* pHamilt,
                     ModulePW::PW_Basis_K* wfc_basis);

    double calne(elecstate::ElecState* pes);

    void itermu(const int iter, elecstate::ElecState* pes);

    void orthog(const int &ik, psi::Psi<std::complex<double>>& psi, Stochastic_WF& stowf);

    void checkemm(const int &ik, const int istep, const int iter, Stochastic_WF& stowf);

    void check_precision(const double ref,const double thr, const std::string info);

    ModuleBase::Chebyshev<double>* p_che = nullptr;

    Stochastic_hchi stohchi;
    Sto_Func<double> stofunc;

	double mu0; // chemical potential; unit in Ry
    bool change;
    double targetne;
    double *spolyv = nullptr;

	public:
    
    int * nchip = nullptr;
    bool check = false;
    double th_ne;
    double KS_ne;
    public:
    int method; //different methods 1: slow, less memory  2: fast, more memory
    ModuleBase::ComplexMatrix* chiallorder = nullptr;
    //chiallorder cost too much memories and should be cleaned after scf.
    void cleanchiallorder();
    //cal shchi = \sqrt{f(\hat{H})}|\chi>
    void calHsqrtchi(Stochastic_WF& stowf);
    //cal Pn = \sum_\chi <\chi|Tn(\hat{h})|\chi>
    void calPn(const int& ik, Stochastic_WF& stowf);
    //cal Tnchi = \sum_n C_n*T_n(\hat{h})|\chi>
    void calTnchi_ik(const int& ik, Stochastic_WF& stowf);
    //cal v^T*M*v
    double vTMv(const double *v, const double * M, const int n);
  private:
    K_Vectors* pkv;

};

#endif// Eelectrons_Iter
