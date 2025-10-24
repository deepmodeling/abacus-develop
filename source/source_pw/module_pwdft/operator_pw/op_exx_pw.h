#ifndef OPEXXPW_H
#define OPEXXPW_H

#include "operator_pw.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/macros.h"
#include "source_base/matrix.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_cell/klist.h"
#include "source_lcao/module_ri/conv_coulomb_pot_k.h"
#include "source_psi/psi.h"
#include "source_base/module_container/ATen/kernels/lapack.h"

#include <memory>
#include <utility>
#include <vector>

namespace hamilt
{

template <typename T, typename Device>
class OperatorEXXPW : public OperatorPW<T, Device>
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    OperatorEXXPW(const int* isk_in,
                  const ModulePW::PW_Basis_K* wfcpw_in,
                  const ModulePW::PW_Basis* rhopw_in,
                  K_Vectors* kv_in,
                  const UnitCell* ucell);

    template <typename T_in, typename Device_in = Device>
    explicit OperatorEXXPW(const OperatorEXXPW<T_in, Device_in> *op_exx);

    virtual ~OperatorEXXPW();

    virtual void act(const int nbands,
                     const int nbasis,
                     const int npol,
                     const T *tmpsi_in,
                     T *tmhpsi,
                     const int ngk_ik = 0,
                     const bool is_first_node = false) const override;

    double cal_exx_energy(psi::Psi<T, Device> *psi_) const;

    void set_psi(psi::Psi<T, Device> &psi_in) const { psi = psi_in; }

    void set_wg(const ModuleBase::matrix *wg_in) { wg = wg_in; }

    void construct_ace() const;

    bool first_iter = true;

    static std::vector<Real> fock_div, erfc_div;

  private:
    const int* isk = nullptr;
    const ModulePW::PW_Basis_K* wfcpw = nullptr;
    const ModulePW::PW_Basis* rhopw = nullptr;
    ModulePW::PW_Basis* rhopw_dev = nullptr; // for device
    const UnitCell *ucell = nullptr;
    Real tpiba = 0;
    
    std::vector<int> get_q_points(const int ik) const;
    const T *get_pw(const int m, const int iq) const;

    void multiply_potential(T *density_recip, int ik, int iq) const;

    void act_op(const int nbands,
                const int nbasis,
                const int npol,
                const T *tmpsi_in,
                T *tmhpsi,
                const int ngk_ik = 0,
                const bool is_first_node = false) const;

    void act_op_ace(const int nbands,
                    const int nbasis,
                    const int npol,
                    const T *tmpsi_in,
                    T *tmhpsi,
                    const int ngk_ik = 0,
                    const bool is_first_node = false) const;

    double cal_exx_energy_op(psi::Psi<T, Device> *psi_) const;

    double cal_exx_energy_ace(psi::Psi<T, Device> *psi_) const;

    void cal_density_recip(const T* psi_nk_real, const T* psi_mq_real, double omega) const;

    void rho_recip2real(const T* rho_recip, T* rho_real, bool add = false, Real factor = 1.0) const;

    mutable int cnt = 0;

    mutable bool potential_got = false;
    
    // pws
//    mutable std::vector<std::unique_ptr<T[]>> pws;

    // k vectors
    K_Vectors *kv = nullptr;

    // psi
    mutable psi::Psi<T, Device> psi;
    const ModuleBase::matrix* wg;

    // real space memory
    T *psi_nk_real = nullptr;
    T *psi_mq_real = nullptr;
    T *density_real = nullptr;
    T *h_psi_real = nullptr;
    // density recip space memory
    T *density_recip = nullptr;
    // h_psi recip space memory
    T *h_psi_recip = nullptr;
    Real *pot = nullptr;

    // Lin Lin's ACE memory, 10.1021/acs.jctc.6b00092
    mutable T* h_psi_ace = nullptr; // H \Psi, W in the paper
    mutable T* psi_h_psi_ace = nullptr; // \Psi^{\dagger} H \Psi, M in the paper
    mutable T* L_ace = nullptr; // cholesky(-M).L, L in the paper
    mutable std::vector<T*> Xi_ace_k; // L^{-1} (H \Psi)^{\dagger}, \Xi in the paper
//    mutable T* Xi_ace = nullptr; // L^{-1} (H \Psi)^{\dagger}, \Xi in the paper

    mutable std::map<int, std::vector<int>> q_points;

    // occupational number
    const ModuleBase::matrix *p_wg;

//    mutable bool update_psi = false;

    Device *ctx = {};
    base_device::DEVICE_CPU* cpu_ctx = {};
    base_device::AbacusDevice_t device = {};

    using ct_Device = typename ct::PsiToContainer<Device>::type;
    using setmem_complex_op = base_device::memory::set_memory_op<T, Device>;
    using setmem_real_op = base_device::memory::set_memory_op<Real, Device>;
    using setmem_real_cpu_op = base_device::memory::set_memory_op<Real, base_device::DEVICE_CPU>;
    using resmem_complex_op = base_device::memory::resize_memory_op<T, Device>;
    using delmem_complex_op = base_device::memory::delete_memory_op<T, Device>;
    using syncmem_complex_op = base_device::memory::synchronize_memory_op<T, Device, Device>;
    using resmem_real_op = base_device::memory::resize_memory_op<Real, Device>;
    using delmem_real_op = base_device::memory::delete_memory_op<Real, Device>;
    using gemm_complex_op = ModuleBase::gemm_op<T, Device>;
    using axpy_complex_op = ModuleBase::axpy_op<T, Device>;
    using vec_add_vec_complex_op = ModuleBase::vector_add_vector_op<T, Device>;
    using dot_op = ModuleBase::dot_real_op<T, Device>;
    using syncmem_complex_c2d_op = base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>;
    using syncmem_complex_d2c_op = base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, Device>;
    using syncmem_real_c2d_op = base_device::memory::synchronize_memory_op<Real, Device, base_device::DEVICE_CPU>;
    using syncmem_real_d2c_op = base_device::memory::synchronize_memory_op<Real, base_device::DEVICE_CPU, Device>;
    using lapack_potrf = container::kernels::lapack_potrf<T, ct_Device>;
    using lapack_trtri = container::kernels::lapack_trtri<T, ct_Device>;

    bool gamma_extrapolation = true;

};

template <typename Real, typename Device>
void get_exx_potential(const K_Vectors* kv,
                       const ModulePW::PW_Basis_K* wfcpw,
                       ModulePW::PW_Basis* rhopw_dev,
                       Real* pot,
                       double tpiba,
                       bool gamma_extrapolation,
                       double ucell_omega,
                       int ik,
                       int iq,
                       bool is_stress = false);

template <typename Real, typename Device>
void get_exx_stress_potential(const K_Vectors* kv,
                              const ModulePW::PW_Basis_K* wfcpw,
                              ModulePW::PW_Basis* rhopw_dev,
                              Real* pot,
                              double tpiba,
                              bool gamma_extrapolation,
                              double ucell_omega,
                              int ik,
                              int iq);

double exx_divergence(Conv_Coulomb_Pot_K::Coulomb_Type coulomb_type,
                      double erfc_omega,
                      const K_Vectors* kv,
                      const ModulePW::PW_Basis_K* wfcpw,
                      ModulePW::PW_Basis* rhopw_dev,
                      double tpiba,
                      bool gamma_extrapolation,
                      double ucell_omega);

} // namespace hamilt

#endif // OPEXXPW_H