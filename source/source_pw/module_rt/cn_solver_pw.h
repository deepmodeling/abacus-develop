#ifndef HSOLVER_CN_SOLVER_PW_H
#define HSOLVER_CN_SOLVER_PW_H

#include "source_base/matrix.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_hsolver/hs_operator.h"
#include "source_hsolver/module_linalg/cn_subspace.h"
#include "source_hsolver/module_linalg/hsolver_linear.h"
#include "source_hsolver/module_linalg/linear_low_rank.h"
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

/** @brief Crank-Nicolson linear solve and preconditioner controls. */
struct CNControls
{
    LinearControls linear;
    PWPreconditioner preconditioner = PWPreconditioner::kinetic;
    bool cn_init = false;
    bool kinetic_enabled = true;
};

LinearMethod parse_linear_method(const std::string& name);
PWPreconditioner parse_pw_precond(const std::string& name);

/** @brief PW Crank-Nicolson numerical solves with fixed previous-step orbitals. */
template <typename T, typename Device>
class CNSolverPW
{
  public:
    /** @brief Configure propagation with explicit numerical controls and orbital layout. */
    CNSolverPW(const ModulePW::PW_Basis_K& basis,
               const CNControls& cn_ctrl,
               int npol,
               const std::vector<int>& global_k,
               const diag_comm_info& comm,
               std::ostream& log);

    /** @brief Invalidate basis-dependent history and restart independent residual checks without releasing buffers.
     *  Call on every rank in the pool after changing the basis or its distribution, including
     *  G-vector ordering or k-point changes that preserve array sizes. Ordinary Hamiltonian
     *  updates and ionic motion at fixed basis do not require this notification.
     */
    void invalidate_basis();
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
               const bool detailed_output);

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
    struct BatchLayout
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
    CNControls cn_ctrl_;
    const int npol_;
    const std::vector<int> global_k_;
    ct::Tensor clean_;
    LinearAlgebra<T, Device> algebra_;
    std::ostream& log_;
    CNSubspace<T, Device> projection_;
    std::vector<KPointState> states_;
    SequenceState sequence_;
    std::unique_ptr<HSolverLinear<T, Device>> linear_solver_;
    ct::Tensor rhs_;
    ct::Tensor inverse_kinetic_;
    ct::Tensor response_workspace_;
    ct::Tensor correction_workspace_;

    void initialize();
    bool tracks_state() const;
    void prepare_sequence(int nk, int ld, int bands, double dt, int step, int iteration);
    bool require_audit(int step, KPointState* state) const;
    SolveDetails solve_kpoint(const LinearOperator<T, Device>& op,
                              const T* previous,
                              T* current,
                              const BatchLayout& batch,
                              int step,
                              int iteration,
                              KPointState* state);
    void retry_kinetic(const LinearOperator<T, Device>& op, T* current, const BatchLayout& batch, LinearSolveResult* result);
    void update_state(KPointState* state, const SolveDetails& details, int step, const BatchLayout& batch, const T* current);
    void report_solve(const SolveDetails& details, int ik, int step, int iteration) const;
    void prepare_buffers(const int nbands, const int nbasis);
    void update_precond(const int ik, const int dim, const T coefficient, const ModuleBase::Vector3<double>& momentum_shift);
};

} // namespace hsolver
#endif
