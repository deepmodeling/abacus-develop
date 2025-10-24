#ifndef DIAGOITERASSIST_H
#define DIAGOITERASSIST_H

#include "source_base/complexmatrix.h"
#include "source_base/macros.h"
#include "source_hamilt/hamilt.h"
#include "source_psi/psi.h"

namespace hsolver
{

template <typename T, typename Device = base_device::DEVICE_CPU>
class DiagoIterAssist
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    static Real PW_DIAG_THR;
    static int PW_DIAG_NMAX;

    static Real LCAO_DIAG_THR;
    static int LCAO_DIAG_NMAX;

    /// average steps of last cg diagonalization for each band.
    static Real avg_iter;
    static bool need_subspace;

    static int SCF_ITER;

    // for psi::Psi structure
     /**
     * @brief Diagonalizes the Hamiltonian in a subspace defined by the given wavefunction.
     *
     * This static function computes the eigenvalues and eigenvectors of the Hamiltonian
     * within the subspace spanned by the provided wavefunction `psi`. The resulting eigenvectors
     * are stored in `evc`, and the corresponding eigenvalues are written to `en`.
     *
     * @tparam T      Data type for computation (e.g., float, double).
     * @tparam Device Device type for computation (e.g., CPU, GPU).
     * @param pHamilt Pointer to the Hamiltonian object.
     * @param psi     Input wavefunction defining the subspace.
     * @param evc     Output container for computed eigenvectors.
     * @param en      Output array for computed eigenvalues.
     * @param n_band  Number of bands (eigenvalues/eigenvectors) to compute. Default is 0 (all).
     * @param is_S_orthogonal If true, assumes the input wavefunction is already orthogonalized.
     */
    static void diag_subspace(const hamilt::Hamilt<T, Device>* const pHamilt,
                              const psi::Psi<T, Device>& psi,
                              psi::Psi<T, Device>& evc,
                              Real *en,
                              int n_band = 0,
                              const bool is_S_orthogonal = false);

    /// @brief use LAPACK to diagonalize the Hamiltonian matrix
    /// @param pHamilt interface to hamiltonian
    /// @param psi wavefunction to diagonalize
    /// @param psi_nr number of rows (nbands)
    /// @param psi_nc number of columns (nbasis)
    /// @param evc new wavefunction
    /// @param en eigenenergies
    /// @note exception handle: if there is no operator initialized in Hamilt, will directly copy value from psi to evc, 
    /// and return all - zero eigenenergies.
    static void diag_subspace_init(
            hamilt::Hamilt<T, Device>* pHamilt,
            const T* psi,
            int psi_nr,
            int psi_nc,
            psi::Psi<T, Device> &evc,
            Real* en,
            const std::function<void(T*, const int)>& add_to_hcc = [](T* null, const int n) {},
            const std::function<void(const T* const, const int, const int)>& export_vcc = [](const T* null, const int n, const int m) {});

    static void diag_heevx(const int nstart,
                            const int nbands,
                            const T *hcc,
                            const int ldh,
                            Real *e,
                            T *vcc);
    static void diag_hegvd(const int nstart,
                            const int nbands,
                            const T *hcc,
                            const T *sc,
                            const int ldh, // nstart
                            Real *e,
                            T *vcc);

    /// @brief calculate Hamiltonian and overlap matrix in subspace spanned by nstart states psi
    /// @param pHamilt : hamiltonian operator carrier
    /// @param psi : wavefunction
    /// @param hcc : Hamiltonian matrix
    /// @param scc : overlap matrix
    static void cal_hs_subspace(const hamilt::Hamilt<T, Device>* pHamilt, // hamiltonian operator carrier
                                                const psi::Psi<T, Device>& psi,     // [in] wavefunction
                                                T *hcc, 
                                                T *scc);

    /// @brief calculate the response matrix from rotation matrix solved by diagonalization of H and S matrix
    /// @param hcc : Hamiltonian matrix
    /// @param scc : overlap matrix
    /// @param nbands : number of bands
    /// @param mat_in : input matrix to be rotated
    /// @param mat_out : output matrix to be rotated
    /// @param mat_col : number of columns of target matrix
    /// @param en : eigenvalues
    static void diag_responce(const T* hcc,
                              const T* scc,
                              const int nbands,
                              const T* mat_in, 
                              T* mat_out, 
                              int mat_col, 
                              Real* en);
    
    /// @brief calculate the response wavefunction psi from rotation matrix solved by diagonalization of H and S matrix
    static void diag_subspace_psi(const T* hcc,
                              const T* scc,
                              const int dim_subspace,
                              psi::Psi<T, Device>& evc,
                              Real* en);

    static bool test_exit_cond(const int& ntry, const int& notconv);

  private:
    constexpr static const Device* ctx = {};

    using hpsi_info = typename hamilt::Operator<T, Device>::hpsi_info;

    using setmem_var_op = base_device::memory::set_memory_op<Real, Device>;
    using resmem_var_op = base_device::memory::resize_memory_op<Real, Device>;
    using delmem_var_op = base_device::memory::delete_memory_op<Real, Device>;
    using syncmem_var_op = base_device::memory::synchronize_memory_op<Real, Device, Device>;
    using syncmem_var_h2d_op
        = base_device::memory::synchronize_memory_op<Real, base_device::DEVICE_GPU, base_device::DEVICE_CPU>;
    using syncmem_var_d2h_op
        = base_device::memory::synchronize_memory_op<Real, base_device::DEVICE_CPU, base_device::DEVICE_GPU>;

    using setmem_complex_op = base_device::memory::set_memory_op<T, Device>;
    using resmem_complex_op = base_device::memory::resize_memory_op<T, Device>;
    using delmem_complex_op = base_device::memory::delete_memory_op<T, Device>;
    using syncmem_complex_op = base_device::memory::synchronize_memory_op<T, Device, Device>;
    using syncmem_complex_h2d_op = base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>;
    using syncmem_complex_d2h_op = base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, Device>;

    static T one;
    static T zero;
};

template <typename T, typename Device>
typename DiagoIterAssist<T, Device>::Real DiagoIterAssist<T, Device>::avg_iter = 0.0;

template <typename T, typename Device>
int DiagoIterAssist<T, Device>::PW_DIAG_NMAX = 30;

template <typename T, typename Device>
typename DiagoIterAssist<T, Device>::Real DiagoIterAssist<T, Device>::PW_DIAG_THR = 1.0e-2;

template <typename T, typename Device>
int DiagoIterAssist<T, Device>::LCAO_DIAG_NMAX = 50;

template <typename T, typename Device>
typename DiagoIterAssist<T, Device>::Real DiagoIterAssist<T, Device>::LCAO_DIAG_THR = 1.0e-12;

template <typename T, typename Device>
bool DiagoIterAssist<T, Device>::need_subspace = false;

template <typename T, typename Device>
int DiagoIterAssist<T, Device>::SCF_ITER = 0;

template <typename T, typename Device>
T DiagoIterAssist<T, Device>::one = static_cast<T>(1.0);

template <typename T, typename Device>
T DiagoIterAssist<T, Device>::zero = static_cast<T>(0.0);
} // namespace hsolver

#endif