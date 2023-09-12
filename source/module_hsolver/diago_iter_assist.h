#ifndef DIAGOITERASSIST_H
#define DIAGOITERASSIST_H

#include "module_base/complexmatrix.h"
#include "module_psi/psi.h"
#include "module_hamilt_general/hamilt.h"

namespace hsolver
{

template<typename FPTYPE, typename Device = psi::DEVICE_CPU>
class DiagoIterAssist
{
  public:
    static FPTYPE PW_DIAG_THR;
    static int PW_DIAG_NMAX;

    /// average steps of last cg diagonalization for each band.
    static FPTYPE avg_iter;
    static bool need_subspace;

    static int SCF_ITER;

    // for CG diagonalization only
    static void diagH_subspace(
        hamilt::Hamilt<FPTYPE, Device>* pHamilt,
        const psi::Psi<std::complex<FPTYPE>, Device> &psi,
        psi::Psi<std::complex<FPTYPE>, Device> &evc,
        FPTYPE *en,
        int n_band = 0);

    static void diagH_subspace_init(
        hamilt::Hamilt<FPTYPE, Device>* pHamilt,
        const std::complex<FPTYPE>* psi,
        int psi_nr,
        int psi_nc,
        psi::Psi<std::complex<FPTYPE>, Device> &evc,
        FPTYPE *en);

    static void diagH_LAPACK(
        const int nstart,
        const int nbands,
        const std::complex<FPTYPE>* hcc,
        const std::complex<FPTYPE>* sc,
        const int ldh, // nstart
        FPTYPE *e,
        std::complex<FPTYPE>* vcc);

    static bool test_exit_cond(const int &ntry, const int &notconv);

  private:
    constexpr static const Device * ctx = {};

    using hpsi_info = typename hamilt::Operator<std::complex<FPTYPE>, Device>::hpsi_info;

    using setmem_var_op = psi::memory::set_memory_op<FPTYPE, Device>;
    using resmem_var_op = psi::memory::resize_memory_op<FPTYPE, Device>;
    using delmem_var_op = psi::memory::delete_memory_op<FPTYPE, Device>;
    using syncmem_var_op = psi::memory::synchronize_memory_op<FPTYPE, Device, Device>;
    using syncmem_var_h2d_op = psi::memory::synchronize_memory_op<FPTYPE, psi::DEVICE_GPU, psi::DEVICE_CPU>;
    using syncmem_var_d2h_op = psi::memory::synchronize_memory_op<FPTYPE, psi::DEVICE_CPU, psi::DEVICE_GPU>;

    using setmem_complex_op = psi::memory::set_memory_op<std::complex<FPTYPE>, Device>;
    using resmem_complex_op = psi::memory::resize_memory_op<std::complex<FPTYPE>, Device>;
    using delmem_complex_op = psi::memory::delete_memory_op<std::complex<FPTYPE>, Device>;
    using syncmem_complex_op = psi::memory::synchronize_memory_op<std::complex<FPTYPE>, Device, Device>;
    using syncmem_complex_h2d_op = psi::memory::synchronize_memory_op<std::complex<FPTYPE>, Device, psi::DEVICE_CPU>;
    using syncmem_complex_d2h_op = psi::memory::synchronize_memory_op<std::complex<FPTYPE>, psi::DEVICE_CPU, Device>;

    static std::complex<FPTYPE> one;
    static std::complex<FPTYPE> zero;
};

template<typename FPTYPE, typename Device>
FPTYPE DiagoIterAssist<FPTYPE, Device>::avg_iter = 0.0;

template<typename FPTYPE, typename Device>
int DiagoIterAssist<FPTYPE, Device>::PW_DIAG_NMAX = 30;

template<typename FPTYPE, typename Device>
FPTYPE DiagoIterAssist<FPTYPE, Device>::PW_DIAG_THR = 1.0e-2;

template<typename FPTYPE, typename Device>
bool DiagoIterAssist<FPTYPE, Device>::need_subspace = false;

template<typename FPTYPE, typename Device>
int DiagoIterAssist<FPTYPE, Device>::SCF_ITER = 0;

template<typename FPTYPE, typename Device>
std::complex<FPTYPE> DiagoIterAssist<FPTYPE, Device>::one = std::complex<FPTYPE>(1, 0);

template<typename FPTYPE, typename Device>
std::complex<FPTYPE> DiagoIterAssist<FPTYPE, Device>::zero = std::complex<FPTYPE>(0, 0);
} // namespace hsolver

#endif