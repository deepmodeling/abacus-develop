#ifndef DIAGOITERASSIST_H
#define DIAGOITERASSIST_H

#include "source_base/macros.h"
#include "source_base/module_device/memory_op.h"
#include "source_hsolver/hs_operator.h"
#include "source_psi/psi.h"

#include <string>

namespace hsolver
{

struct diag_comm_info;

template <typename T, typename Device = base_device::DEVICE_CPU>
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

    /**
     * @brief Diagonalizes H in the subspace spanned by nstart vectors.
     *
     * Builds the nstart*nstart matrices <psi|H|psi> (and <psi|S|psi> unless the
     * input is S-orthogonal), solves the small eigenproblem and rotates psi
     * into the lowest n_band eigenvectors, written to evc.
     *
     * @param op      applies H and S to block vectors
     * @param psi     [in]  nstart vectors, leading dimension dmax
     * @param evc     [out] n_band vectors, leading dimension dmax; may alias psi
     * @param nstart  number of input vectors
     * @param n_band  number of eigenvectors wanted (<= nstart)
     * @param dmin    active length of each vector
     * @param dmax    leading dimension of psi and evc
     * @param en      [out] n_band eigenvalues (host memory)
     * @param is_S_orthogonal if true, psi is already S-orthonormal and the
     *        standard eigenproblem is solved instead of the generalized one
     */
    static void diag_subspace(const HSOperator<T, Device>& op,
                              const T* psi,
                              T* evc,
                              const int nstart,
                              const int n_band,
                              const int dmin,
                              const int dmax,
                              Real* en,
                              const diag_comm_info& diag_comm,
                              const bool is_S_orthogonal = false);

    /// psi::Psi flavour of diag_subspace(): nstart = psi.get_nbands(),
    /// n_band = 0 means all of them, dimensions taken from psi.
    static void diag_subspace(const HSOperator<T, Device>& op,
                              const psi::Psi<T, Device>& psi,
                              psi::Psi<T, Device>& evc,
                              Real* en,
                              const diag_comm_info& diag_comm,
                              int n_band = 0,
                              const bool is_S_orthogonal = false);

    /// @brief subspace diagonalization used to build the starting wavefunction
    /// @param op interface to H and S; op.add_to_subspace_h() and
    /// op.export_subspace_vec() are called around the small eigenproblem
    /// @param psi vectors spanning the subspace
    /// @param psi_nr number of rows (nbands)
    /// @param psi_nc number of columns (nbasis)
    /// @param evc new wavefunction
    /// @param en eigenenergies
    /// @param basis_type "lcao", "lcao_in_pw" or "pw"; together with calculation it selects
    /// how the rotation matrix is applied to psi
    /// @param calculation "scf", "nscf", "md", "relax", ...
    static void diag_subspace_init(const HSOperator<T, Device>& op,
                                   const T* psi,
                                   int psi_nr,
                                   int psi_nc,
                                   psi::Psi<T, Device>& evc,
                                   Real* en,
                                   const std::string& basis_type,
                                   const std::string& calculation,
                                   const diag_comm_info& diag_comm);

    static void diag_heevx(const int nstart,
                            const int nbands,
                            const T *hcc,
                            const int ldh,
                            Real *e,
                            T *vcc);
    static void diag_hegvd(const int nstart,
                            const int nbands,
                            const T *hcc,
                            T *sc,
                            const int ldh, // nstart
                            Real *e,
                            T *vcc);

    /// @brief calculate Hamiltonian and overlap matrix in subspace spanned by nstart states psi
    /// @param op : applies H and S
    /// @param psi : wavefunction
    /// @param hcc : Hamiltonian matrix
    /// @param scc : overlap matrix
    static void cal_hs_subspace(const HSOperator<T, Device>& op,
                                const psi::Psi<T, Device>& psi, // [in] wavefunction
                                T* hcc,
                                T* scc,
                                const diag_comm_info& diag_comm);

    /// @brief calculate the response matrix from rotation matrix solved by diagonalization of H and S matrix
    /// @param hcc : Hamiltonian matrix
    /// @param scc : overlap matrix
    /// @param nbands : number of bands
    /// @param mat_in : input matrix to be rotated
    /// @param mat_out : output matrix to be rotated
    /// @param mat_col : number of columns of target matrix
    /// @param en : eigenvalues
    static void diag_responce(const T* hcc,
                              T* scc,
                              const int nbands,
                              const T* mat_in, 
                              T* mat_out, 
                              int mat_col, 
                              Real* en);
    
    /// @brief calculate the response wavefunction psi from rotation matrix solved by diagonalization of H and S matrix
    static void diag_subspace_psi(const T* hcc,
                              T* scc,
                              const int dim_subspace,
                              psi::Psi<T, Device>& evc,
                              Real* en);

  private:
    constexpr static const Device* ctx = {};

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
bool DiagoIterAssist<T, Device>::need_subspace = false;

template <typename T, typename Device>
int DiagoIterAssist<T, Device>::SCF_ITER = 0;

template <typename T, typename Device>
T DiagoIterAssist<T, Device>::one = static_cast<T>(1.0);

template <typename T, typename Device>
T DiagoIterAssist<T, Device>::zero = static_cast<T>(0.0);
} // namespace hsolver

#endif
