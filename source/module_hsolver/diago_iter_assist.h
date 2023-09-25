#ifndef DIAGOITERASSIST_H
#define DIAGOITERASSIST_H

#include "module_base/complexmatrix.h"
#include "module_psi/psi.h"
#include "module_hamilt_general/hamilt.h"

#include <module_base/macros.h>
namespace hsolver
{

template<typename T, typename Device = psi::DEVICE_CPU>
class DiagoIterAssist
{
  private:
    using Real = typename GetTypeReal<T>::type;
  public:
    static Real PW_DIAG_THR;
    static int PW_DIAG_NMAX;

    /// average steps of last cg diagonalization for each band.
    static Real avg_iter;
    static bool need_subspace;

    static int SCF_ITER;

    // for CG diagonalization only
    static void diagH_subspace(
        hamilt::Hamilt<T, Device>* pHamilt,
        const psi::Psi<T, Device> &psi,
        psi::Psi<T, Device> &evc,
        Real *en,
        int n_band = 0);

    static void diagH_subspace_init(
        hamilt::Hamilt<T, Device>* pHamilt,
        const T* psi,
        int psi_nr,
        int psi_nc,
        psi::Psi<T, Device> &evc,
        Real *en);

    static void diagH_LAPACK(
        const int nstart,
        const int nbands,
        const T* hcc,
        const T* sc,
        const int ldh, // nstart
        Real *e,
        T* vcc);

    static bool test_exit_cond(const int &ntry, const int &notconv);

  private:
    constexpr static const Device * ctx = {};

    using hpsi_info = typename hamilt::Operator<T, Device>::hpsi_info;

    using setmem_var_op = psi::memory::set_memory_op<Real, Device>;
    using resmem_var_op = psi::memory::resize_memory_op<Real, Device>;
    using delmem_var_op = psi::memory::delete_memory_op<Real, Device>;
    using syncmem_var_op = psi::memory::synchronize_memory_op<Real, Device, Device>;
    using syncmem_var_h2d_op = psi::memory::synchronize_memory_op<Real, psi::DEVICE_GPU, psi::DEVICE_CPU>;
    using syncmem_var_d2h_op = psi::memory::synchronize_memory_op<Real, psi::DEVICE_CPU, psi::DEVICE_GPU>;

    using setmem_complex_op = psi::memory::set_memory_op<T, Device>;
    using resmem_complex_op = psi::memory::resize_memory_op<T, Device>;
    using delmem_complex_op = psi::memory::delete_memory_op<T, Device>;
    using syncmem_complex_op = psi::memory::synchronize_memory_op<T, Device, Device>;
    using syncmem_complex_h2d_op = psi::memory::synchronize_memory_op<T, Device, psi::DEVICE_CPU>;
    using syncmem_complex_d2h_op = psi::memory::synchronize_memory_op<T, psi::DEVICE_CPU, Device>;

    static T one;
    static T zero;
};

template<typename T, typename Device>
typename DiagoIterAssist<T, Device>::Real DiagoIterAssist<T, Device>::avg_iter = 0.0;

template<typename T, typename Device>
int DiagoIterAssist<T, Device>::PW_DIAG_NMAX = 30;

template<typename T, typename Device>
typename DiagoIterAssist<T, Device>::Real DiagoIterAssist<T, Device>::PW_DIAG_THR = 1.0e-2;

template<typename T, typename Device>
bool DiagoIterAssist<T, Device>::need_subspace = false;

template<typename T, typename Device>
int DiagoIterAssist<T, Device>::SCF_ITER = 0;

template<typename T, typename Device>
T DiagoIterAssist<T, Device>::one = static_cast<T>(1.0);

template<typename T, typename Device>
T DiagoIterAssist<T, Device>::zero = static_cast<T>(0.0);
} // namespace hsolver

#endif