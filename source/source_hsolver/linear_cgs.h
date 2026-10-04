#ifndef HSOLVER_LINEAR_CGS_H
#define HSOLVER_LINEAR_CGS_H
#include "source_hsolver/linear_workspace.h"

namespace hsolver
{
/** @brief Right-preconditioned CGS with a compact active set of columns. */
template <typename T, typename Device = base_device::DEVICE_CPU>
class LinearCGS final
{
  private:
    enum WorkspaceSlot
    {
        residual_slot = 0,
        shadow_slot = 1,
        direction_slot = 2,
        q_slot = 3,
        u_slot = 4,
        image_slot = 5,
        scratch_slot = 6,
        update_slot = 7,
        solution_slot = 8
    };
    const double tolerance_;
    const int max_iter_;
    LinearWorkspace<T, Device> work_;
    int ld_ = 0;
    int dim_ = 0;
    int active_ = 0;
    std::vector<T> beta_;
    std::vector<T> rho_;
    std::vector<T> rho_prev_;
    std::vector<T> alpha_;
    std::vector<T> denominator_;
    std::vector<T> norm_;
    std::vector<double> threshold_;
    std::vector<int> original_;
    std::vector<int> swaps_;

  public:
    LinearCGS(const double tolerance, const int max_iter, const diag_comm_info& comm);
    /**
     * @brief Solve A*x=b in Device memory, preserving padding and caller column order.
     * @note Uses tolerance*norm(b), or tolerance for a zero right-hand side.
     */
    LinearSolveResult solve(const LinearOperator<T, Device>& op,
                            const LinearOperator<T, Device>& preconditioner,
                            const int ld,
                            const int nband,
                            const int dim,
                            T* x,
                            const T* b);

  private:
    bool iterate(const LinearOperator<T, Device>& op,
                 const LinearOperator<T, Device>& preconditioner,
                 const bool first_iteration,
                 LinearSolveResult* result);
    void retire_converged();
};
} // namespace hsolver
#endif
