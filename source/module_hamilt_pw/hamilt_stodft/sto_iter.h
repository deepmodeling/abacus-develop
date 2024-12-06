#ifndef STO_ITER_H
#define STO_ITER_H
#include "module_base/math_chebyshev.h"
#include "module_elecstate/elecstate_pw.h"
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
  private:
    using Real = typename GetTypeReal<T>::type;
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
              StoChe<Real, Device>& stoche,
              hamilt::HamiltSdftPW<T, Device>* p_hamilt_sto);

    /**
     * @brief sum demet and eband energies for each k point and each band
     * 
     * @param stowf stochastic wave function
     * @param pes elecstate
     * @param pHamilt hamiltonian
     * @param wfc_basis wfc pw basis
     */
    void sum_stoeband(Stochastic_WF<T, Device>& stowf,
                     elecstate::ElecStatePW<T, Device>* pes,
                     hamilt::Hamilt<T, Device>* pHamilt,
                     ModulePW::PW_Basis_K* wfc_basis);

    /**
     * @brief calculate the density
     * 
     * @param ucell reference to unit cell
     * @param stowf stochastic wave function
     * @param pes elecstate
     * @param wfc_basis wfc pw basis
     */
    void cal_storho(const UnitCell& ucell,
                    Stochastic_WF<T, Device>& stowf,
                    elecstate::ElecStatePW<T, Device>* pes,
                    ModulePW::PW_Basis_K* wfc_basis);

    /**
     * @brief calculate total number of electrons
     * 
     * @param pes elecstate
     * @return double 
     */
    double calne(elecstate::ElecState* pes);

    /**
     * @brief solve ne(mu) = ne_target and get chemical potential mu
     * 
     * @param iter scf iteration index
     * @param pes elecstate
     */
    void itermu(const int iter, elecstate::ElecState* pes);

    /**
     * @brief orthogonalize stochastic wave functions with KS wave functions
     * 
     * @param ik k point index
     * @param psi KS wave functions
     * @param stowf stochastic wave functions
     */
    void orthog(const int& ik, psi::Psi<T, Device>& psi, Stochastic_WF<T, Device>& stowf);

    /**
     * @brief check emax and emin
     * 
     * @param ik k point index
     * @param istep ion step index
     * @param iter scf iteration index
     * @param stowf stochastic wave functions
     */
    void checkemm(const int& ik, const int istep, const int iter, Stochastic_WF<T, Device>& stowf);

    /**
     * @brief check precision of Chebyshev expansion
     * 
     * @param ref reference value
     * @param thr threshold
     * @param info information
     */
    void check_precision(const double ref, const double thr, const std::string info);

    ModuleBase::Chebyshev<double, Device>* p_che = nullptr;

    Sto_Func<double> stofunc;
    hamilt::HamiltSdftPW<T, Device>* p_hamilt_sto = nullptr;

    double mu0; // chemical potential; unit in Ry
    bool change;
    double targetne;
    Real* spolyv = nullptr;     //[Device] coefficients of Chebyshev expansion
    Real* spolyv_cpu = nullptr; //[CPU] coefficients of Chebyshev expansion

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
    /**
     * @brief return cpu dot result
     * @param x [Device]
     * @param y [Device]
     * @param result [CPU] dot result
     */
    void dot(const int& n, const Real* x, const int& incx, const Real* y, const int& incy,  Real& result);
  private:
    const Device* ctx = {};
    const base_device::DEVICE_CPU* cpu_ctx = {};
    using ct_Device = typename container::PsiToContainer<Device>::type;
    using setmem_var_op = base_device::memory::set_memory_op<Real, Device>;
    using syncmem_var_h2d_op = base_device::memory::synchronize_memory_op<Real, Device, base_device::DEVICE_CPU>;
    using syncmem_var_d2h_op = base_device::memory::synchronize_memory_op<Real, base_device::DEVICE_CPU, Device>;
    using cpymem_complex_op = base_device::memory::synchronize_memory_op<T, Device, Device>;
    using resmem_var_op = base_device::memory::resize_memory_op<Real, Device>;
    using delmem_var_op = base_device::memory::delete_memory_op<Real, Device>;
    using resmem_complex_op = base_device::memory::resize_memory_op<T, Device>;
    using delmem_complex_op = base_device::memory::delete_memory_op<T, Device>;
    using castmem_d2z_op = base_device::memory::cast_memory_op<T, Real, Device, Device>;
    using castmem_var_d2h_op = base_device::memory::cast_memory_op<double, Real, base_device::DEVICE_CPU, Device>;
    using gemv_op = hsolver::gemv_op<T, Device>;
};

#endif // Eelectrons_Iter
