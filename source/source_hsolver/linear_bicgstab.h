#ifndef MODULE_HSOLVER_LINEAR_BICGSTAB_H_
#define MODULE_HSOLVER_LINEAR_BICGSTAB_H_
#include "source_hsolver/linear_workspace.h"
namespace hsolver
{
/** @brief Right-preconditioned BiCGSTAB for independent complex right-hand sides. */
template <typename T, typename Device = base_device::DEVICE_CPU>
class LinearBiCGSTAB final
{
  private:
    enum WorkspaceSlot
    {
        residual_slot = 0,
        shadow_slot = 1,
        direction_slot = 2,
        image_slot = 3,
        intermediate_slot = 4,
        precond_direction_slot = 5,
        precond_intermediate_slot = 6,
        intermediate_image_slot = 7,
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
    std::vector<T> omega_;
    std::vector<T> denominator_;
    std::vector<T> tt_;
    std::vector<T> norm_;
    std::vector<double> threshold_;
    std::vector<int> original_;
    std::vector<int> swaps_;
    std::vector<int> skip_;

  public:
    LinearBiCGSTAB(const double tolerance, const int max_iter, const diag_comm_info& comm);
    /**
     * @brief Solve A*x=b, using x as initial guess and M as inverse preconditioner.
     * @param ld Column stride, at least dim. Padding is neither read nor overwritten.
     * @param dim Local valid row count; may be zero on an MPI participant.
     * @note Accepts only true residuals below tolerance*max(1, norm(b)) per column.
     */
    LinearSolveResult solve(const LinearOperator<T, Device>& op,
                            const LinearOperator<T, Device>& preconditioner,
                            const int ld,
                            const int nband,
                            const int dim,
                            T* x,
                            const T* b);

  private:
    bool iterate(const LinearOperator<T, Device>& op, const LinearOperator<T, Device>& preconditioner, LinearSolveResult* result);
    void retire(const int band);
    void retire_converged(const bool alpha_step, const bool identity);
    bool breakdown(const T divisor, const int band, LinearSolveResult* result) const;
};
} // namespace hsolver
#endif
