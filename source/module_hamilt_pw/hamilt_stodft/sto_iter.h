#ifndef STO_ITER_H
#define STO_ITER_H
#include "module_base/math_chebyshev.h"
#include "module_elecstate/elecstate.h"
#include "module_hamilt_general/hamilt.h"
#include "module_hamilt_pw/hamilt_stodft/hamilt_sdft_pw.h"
#include "module_psi/psi.h"
#include "sto_che.h"
#include "sto_func.h"
#include "sto_wf.h"

//----------------------------------------------
// Solve for the new electron density and iterate
// until SCF loop converges (mu, epsilon, rho)
// mu: chemical potential
// epsilon: eigenvalues
// rho: charge density
//----------------------------------------------

template <typename T, typename Device = base_device::DEVICE_CPU>
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
     * @param stoche Chebyshev expansion for sDFT
     * @param p_hamilt_sto hamiltonian for sDFT
     *
     */
    void init(K_Vectors* pkv_in,
              ModulePW::PW_Basis_K* wfc_basis,
              Stochastic_WF<T, Device>& stowf,
              StoChe<double>& stoche,
              hamilt::HamiltSdftPW<T, Device>* p_hamilt_sto);

    void sum_stoband(Stochastic_WF<T, Device>& stowf,
                     elecstate::ElecState* pes,
                     hamilt::Hamilt<T, Device>* pHamilt,
                     ModulePW::PW_Basis_K* wfc_basis);

    double calne(elecstate::ElecState* pes);

    void itermu(const int iter, elecstate::ElecState* pes);

    void orthog(const int& ik, psi::Psi<T, Device>& psi, Stochastic_WF<T, Device>& stowf);

    void checkemm(const int& ik, const int istep, const int iter, Stochastic_WF<T, Device>& stowf);

    void check_precision(const double ref, const double thr, const std::string info);

    ModuleBase::Chebyshev<double>* p_che = nullptr;

    Sto_Func<double> stofunc;
    hamilt::HamiltSdftPW<T, Device>* p_hamilt_sto = nullptr;

    double mu0; // chemical potential; unit in Ry
    bool change;
    double targetne;
    double* spolyv = nullptr;

  public:
    int* nchip = nullptr;
    bool check = false;
    double th_ne;
    double KS_ne;

  public:
    int method; // different methods 1: slow, less memory  2: fast, more memory
    // cal shchi = \sqrt{f(\hat{H})}|\chi>
    void calHsqrtchi(Stochastic_WF<T, Device>& stowf);
    // cal Pn = \sum_\chi <\chi|Tn(\hat{h})|\chi>
    void calPn(const int& ik, Stochastic_WF<T, Device>& stowf);
    // cal Tnchi = \sum_n C_n*T_n(\hat{h})|\chi>
    void calTnchi_ik(const int& ik, Stochastic_WF<T, Device>& stowf);

  private:
    K_Vectors* pkv;
};

#endif // Eelectrons_Iter
