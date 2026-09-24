#ifndef OPEXXPW_H
#define OPEXXPW_H

#include "op_pw.h"
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
                  const UnitCell* ucell,
                  const bool separate_loop_in,
                  const Real hybrid_alpha_in,
                  const CoulombParam& coulomb_param_in);

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
    bool separate_loop = false;
    Real hybrid_alpha = 0.0;
    CoulombParam coulomb_param;

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

    void act_op_kpar(const int nbands,
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
    // cache of psi_nk in real space for all bands on the active grid
    mutable T* psi_nk_real_cache = nullptr;
    mutable int psi_nk_cache_size = 0; // number of bands currently allocated in the cache

    // ---- EXX FFT grid + band-batch layer ----
    // All EXX PW entry points (act_op, act_op_kpar,
    // cal_exx_energy_op) share one code path built on these primitives; the
    // batched kernels are a specialization of the per-band operations,
    // selected inside the primitives (exx_grid_active). The active real-space
    // grid is the small ecut_exx box when usable, else the wfcpw box;
    // callers must not branch on this themselves.
    void maybe_setup_exx_grid() const;
    // whether the map-based (batched-kernel) path can run at all: the small
    // grid is usable, or the full box is local to this rank
    bool exx_grid_active() const;
    // band chunk width of one batched round, 1..nbands: the exx_batch_size
    // INPUT (0 = all bands). Only meaningful when exx_grid_active(); the
    // result is independent of the chunking.
    int exx_band_chunk(const int nbands) const;
    int exx_grid_size() const { return exx_sg_ok ? sg_nxyz : wfcpw->nrxx; }
    // G-vector -> active FFT box index maps, hiding the grid choice: the
    // small-grid maps when exx_sg_ok, the basis's own ig2ixyz arrays on CUDA,
    // and the host-built full-box maps (full_map_*) on CPU otherwise
    const int* active_map_rho() const;
    const int* active_map_wfc(const int ik) const;
    // fill psi_nk_real_cache: all bands of psi_in at k-point ik, real space on
    // the active grid
    void cache_psi_nk_real(const int nbands, const int nbasis, const T* psi_in, const int ik) const;
    // psi_mq_real: single band psi_m at k-point iq, real space on the active
    // grid (psi_m points to a band block with row stride psi_stride)
    void wfc_to_real_exx(const T* psi_m, const int iq, const int psi_stride) const;
    // tmhpsi += factor * V_x|psi> for all bands; needs psi_nk_real_cache and
    // psi_mq_real filled, and `pot` holding the Coulomb kernel for (ik, iq)
    void apply_fock_all_bands(const int nbands, const int nbasis, const int iq, const Real factor, T* tmhpsi) const;
    // pair densities psi_nk* psi_mq of bands [start, start+count) in the
    // rhopw_dev G-space; pair_density(n) then returns band n (global index)
    // without recomputation
    void prepare_pair_densities(const int start, const int count) const;
    const T* pair_density(const int n) const;

    // Small EXX FFT grid sized by ecut_exx (like QE's dfftt): when usable, the
    // batched path runs all FFTs on this smaller grid instead of the wfcpw grid.
    // Requires all |k+G|^2 < ecut_exx (otherwise the wavefunctions do not fit).
    // Available on CPU and CUDA; on ROCm the GPU path keeps the full grid
    // (exx_sg_ok stays false there).
    void setup_exx_small_grid() const;
    // host-built full-box maps for the batched path on CPU when the small
    // grid is not active (the GPU bases own these as ig2ixyz*)
    void setup_full_grid_maps() const;
    mutable bool exx_sg_init = false;
    mutable bool exx_sg_ok = false;
    mutable int sg_nx = 0, sg_ny = 0, sg_nz = 0, sg_nxyz = 0;
    mutable int* sg_map_rho = nullptr; // rhopw_dev G-vector ig -> small box index [rhopw_dev->npw]
    mutable int* sg_map_wfc = nullptr; // wfcpw (G+k) ig -> small box index [npwk_max * nks]
    mutable int* full_map_rho = nullptr; // rhopw_dev G-vector ig -> full box index (CPU only)
    mutable int* full_map_wfc = nullptr; // wfcpw (G+k) ig -> full box index (CPU only)

    // Batched n-band loop of act_op / act_op_kpar: all bands of one
    // (iq, m_iband) pair are processed with batched FFTs and a handful of
    // kernel launches instead of ~8 launches per band. Requires the psi_nk
    // real-space cache and a local (non-distributed) FFT box. `factor` is the
    // full accumulation weight (including hybrid_alpha and the k/q weights).
    void apply_exx_nbatched(const int nbands,
                            const int nbasis,
                            const T* psi_mq_real,
                            const Real factor,
                            T* tmhpsi) const;

    // reciprocal-space density of bands [0, nbands) of the nk_real block
    // (rhopw_dev G-space), ends up in dens_pw_batch
    void calc_density_pw_nbatched(const int nbands, const T* nk_real, const T* psi_mq_real) const;

    mutable void* exx_fft_plan = nullptr;   // batched FFT plan (void* to keep FFTW/cuFFT out of the header)
    mutable void* exx_fft_plan1 = nullptr;  // batch-1 FFT plan (psi_mq on the active grid)
    mutable int exx_fft_plan_nx = 0;        // grid dims the plans were created for
    mutable int exx_fft_plan_ny = 0;
    mutable int exx_fft_plan_nz = 0;
    mutable T* dens_box_batch = nullptr;    // band-chunk * nxyz box buffer
    mutable T* dens_pw_batch = nullptr;     // band-chunk * npw (rhopw_dev) plane-wave buffer
    mutable int exx_batch_alloc = 0;        // band chunk width the buffers/plan are allocated for
    mutable int dens_chunk_base = 0;        // global band index of dens_pw_batch[0]

    // create/recreate the batched FFT plans and work buffers for the active grid
    void ensure_exx_batch(const int nbands) const;
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
    using vec_add_vec_complex_op = ModuleBase::vector_add_vector_op<T, Device, Real>;
    using dot_op = ModuleBase::dot_real_op<T, Device>;
    using syncmem_complex_c2d_op = base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>;
    using syncmem_complex_d2c_op = base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, Device>;
    using syncmem_real_c2d_op = base_device::memory::synchronize_memory_op<Real, Device, base_device::DEVICE_CPU>;
    using syncmem_real_d2c_op = base_device::memory::synchronize_memory_op<Real, base_device::DEVICE_CPU, Device>;
    using resmem_int_op = base_device::memory::resize_memory_op<int, Device>;
    using delmem_int_op = base_device::memory::delete_memory_op<int, Device>;
    using syncmem_int_h2d_op = base_device::memory::synchronize_memory_op<int, Device, base_device::DEVICE_CPU>;
    using syncmem_int_d2h_op = base_device::memory::synchronize_memory_op<int, base_device::DEVICE_CPU, Device>;
    using lapack_potrf = container::kernels::lapack_potrf<T, ct_Device>;
    using lapack_trtri = container::kernels::lapack_trtri<T, ct_Device>;

    bool gamma_extrapolation = true;

};

// Explicit specializations must be declared before any implicit instantiation.
// The extern template declarations below would otherwise instantiate the
// generic cal_density_recip / rho_recip2real members.
template <>
void OperatorEXXPW<std::complex<double>, base_device::DEVICE_CPU>::cal_density_recip(
    const std::complex<double>* psi_nk_real,
    const std::complex<double>* psi_mq_real,
    double omega) const;

template <>
void OperatorEXXPW<std::complex<float>, base_device::DEVICE_CPU>::cal_density_recip(
    const std::complex<float>* psi_nk_real,
    const std::complex<float>* psi_mq_real,
    double omega) const;

template <>
void OperatorEXXPW<std::complex<double>, base_device::DEVICE_CPU>::rho_recip2real(
    const std::complex<double>* rho_recip,
    std::complex<double>* rho_real,
    bool add,
    double factor) const;

template <>
void OperatorEXXPW<std::complex<float>, base_device::DEVICE_CPU>::rho_recip2real(
    const std::complex<float>* rho_recip,
    std::complex<float>* rho_real,
    bool add,
    float factor) const;

#if ((defined __CUDA) || (defined __ROCM))
template <>
void OperatorEXXPW<std::complex<double>, base_device::DEVICE_GPU>::cal_density_recip(
    const std::complex<double>* psi_nk_real,
    const std::complex<double>* psi_mq_real,
    double omega) const;

template <>
void OperatorEXXPW<std::complex<float>, base_device::DEVICE_GPU>::cal_density_recip(
    const std::complex<float>* psi_nk_real,
    const std::complex<float>* psi_mq_real,
    double omega) const;

template <>
void OperatorEXXPW<std::complex<double>, base_device::DEVICE_GPU>::rho_recip2real(
    const std::complex<double>* rho_recip,
    std::complex<double>* rho_real,
    bool add,
    double factor) const;

template <>
void OperatorEXXPW<std::complex<float>, base_device::DEVICE_GPU>::rho_recip2real(
    const std::complex<float>* rho_recip,
    std::complex<float>* rho_real,
    bool add,
    float factor) const;
#endif

extern template class OperatorEXXPW<std::complex<float>, base_device::DEVICE_CPU>;
extern template class OperatorEXXPW<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
extern template class OperatorEXXPW<std::complex<float>, base_device::DEVICE_GPU>;
extern template class OperatorEXXPW<std::complex<double>, base_device::DEVICE_GPU>;
#endif

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
                       bool is_stress,
                       const CoulombParam& coulomb_param_in);

template <typename Real, typename Device>
void get_exx_stress_potential(const K_Vectors* kv,
                              const ModulePW::PW_Basis_K* wfcpw,
                              ModulePW::PW_Basis* rhopw_dev,
                              Real* pot,
                              double tpiba,
                              bool gamma_extrapolation,
                              double ucell_omega,
                              int ik,
                              int iq,
                              const CoulombParam& coulomb_param_in);

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
