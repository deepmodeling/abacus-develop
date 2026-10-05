#ifndef HSOLVER_PW_TDDFT_H
#define HSOLVER_PW_TDDFT_H

#include "source_base/matrix.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_hsolver/cn_subspace.h"
#include "source_hsolver/hs_operator.h"
#include "source_hsolver/hsolver_linear.h"
#include "source_hsolver/linear_low_rank.h"
#include "source_hsolver/orthonormal.h"
#include "source_psi/psi.h"

#include <iosfwd>
#include <memory>
#include <string>

namespace hsolver
{

enum class PWPreconditioner
{
    none,
    kinetic,
    kinetic_recycle,
    kinetic_subspace
};

/** @brief Numerical and propagation options, resolved once at the input boundary. */
struct PWLinearOptions
{
    OrthMethod orthonormal = OrthMethod::cholesky;
    LinearSolveOptions linear;
    PWPreconditioner preconditioner = PWPreconditioner::kinetic;
    bool cn_init = false;
    bool kinetic_enabled = true;
};

LinearMethod parse_linear_method(const std::string& name);
/** @brief Step diagnostics, replicated within each k-point pool. */
struct TDOrthStats
{
    double before = 0.0;
    double after = 0.0;
    double orth_seconds = 0.0;
    int calls = 0;
    int passes = 0;
    int fallbacks = 0;
    int skipped = 0;
    int orth_warnings = 0;
    std::array<int, 8> events{};
    std::string reason;
};
PWPreconditioner parse_pw_precond(const std::string& name);

/** @brief PW Crank-Nicolson solves and endpoint Hamiltonian expectations. */
template <typename T, typename Device>
class HSolverPWTDDFT
{
  public:
    /** @brief Configure propagation using typed, explicit solver options. */
    HSolverPWTDDFT(const ModulePW::PW_Basis_K& basis, const PWLinearOptions& options, const diag_comm_info& comm, std::ostream& log);

    /** @brief Invalidate basis-dependent history and restart independent residual checks without releasing buffers.
     *  Call on every rank in the pool after changing the basis or its distribution, including
     *  G-vector ordering or k-point changes that preserve array sizes. Ordinary Hamiltonian
     *  updates and ionic motion at fixed basis do not require this notification.
     */
    void invalidate_basis();
    /** @brief Reset diagnostics without invalidating preconditioner history. */
    void reset_orth_stats();
    /** @brief Correct initial orbitals before rebuilding occupations and density. */
    void orthonormalize(psi::Psi<T, Device>* current);
    /** @brief Occupation-weighted electron count for the local k-point pool. */
    double wave_electrons(const ModuleBase::matrix& occupations) const;
    const TDOrthStats& orth_stats() const
    {
        return orth_stats_;
    }

    /** @brief Propagate from the fixed previous step, retaining the current iterate as the initial guess.
     *  @param dt Electronic time step in Hartree atomic units.
     *  @param momentum_shift Propagation vector potential in Hartree atomic units (inverse Bohr).
     *  @param istep Electronic propagation step, increasing across electronic substeps within an MD step.
     *  @param iter SCF iteration within this electronic step.
     */
    void solve(HSOperator<T, Device>& op,
               const psi::Psi<T, Device>& previous,
               psi::Psi<T, Device>* current,
               const double dt,
               const ModuleBase::Vector3<double>& momentum_shift,
               const int istep,
               const int iter,
               const bool detailed_output,
               std::ostream& log);

    /** @brief Evaluate band expectations after the caller restores the endpoint Hamiltonian. */
    void cal_band_energy(HSOperator<T, Device>& op, const psi::Psi<T, Device>& current, ModuleBase::matrix* energies);

  private:
    using Real = typename GetTypeReal<T>::type;
    struct KPointState
    {
        LinearResponse<T, Device> response;
        unsigned int solve_count = 0;
        int independent_step = -1;
    };
    struct SequenceState
    {
        double dt = -1.0;
        int step = -1;
        int iteration = -1;
        int ld = -1;
        int bands = -1;
    };
    struct SolveBatch
    {
        int ld;
        int dim;
        int bands;
    };
    struct SolveDetails
    {
        LinearSolveResult linear;
        int coarse_rank = 0;
        bool projected = false;
        bool cn_initial = false;
        bool retried = false;
    };
    const ModulePW::PW_Basis_K& basis_;
    const diag_comm_info comm_;
    PWLinearOptions options_;
    LinearAlgebra<T, Device> algebra_;
    Orthonormal<T, Device> orthonormal_;
    TDOrthStats orth_stats_;
    std::vector<std::vector<double>> orth_norms_;
    CNSubspace<T, Device> projection_;
    std::vector<KPointState> states_;
    SequenceState sequence_;
    std::unique_ptr<HSolverLinear<T, Device>> linear_solver_;
    LinearWorkspace<T, Device> band_products_;
    ct::Tensor rhs_;
    ct::Tensor hpsi_;
    ct::Tensor inverse_kinetic_;
    ct::Tensor response_workspace_;
    ct::Tensor correction_workspace_;

    void initialize(std::ostream& log);
    bool tracks_state() const;
    void prepare_sequence(int nk, int ld, int bands, double dt, int step, int iteration);
    bool require_audit(int step, KPointState* state) const;
    SolveDetails solve_kpoint(const LinearOperator<T, Device>& op,
                              const T* previous,
                              T* current,
                              const SolveBatch& batch,
                              int step,
                              int iteration,
                              KPointState* state);
    void retry_kinetic(const LinearOperator<T, Device>& op, T* current, const SolveBatch& batch, LinearSolveResult* result);
    void update_state(KPointState* state, const SolveDetails& details, int step, const SolveBatch& batch, const T* current);
    void report_solve(const SolveDetails& details, int ik, int step, int iteration, double elapsed, std::ostream& log) const;
    void prepare_buffers(const int nbands, const int nbasis);
    void update_precond(const int ik, const int dim, const T coefficient, const ModuleBase::Vector3<double>& momentum_shift);
    void correct_orbitals(T* current, int ld, int dim, int bands, int ik);
};

} // namespace hsolver
#endif
