#ifndef HSOLVER_LINEAR_GMRES_H
#define HSOLVER_LINEAR_GMRES_H
#include "source_hsolver/linear_algebra.h"
#include "source_hsolver/linear_workspace.h"

namespace hsolver
{
/** @brief Restarted right-preconditioned GMRES with independent batched columns. */
template <typename T, typename Device>
class LinearGMRES
{
  private:
    using Wide = std::complex<double>;
    using Real = typename GetTypeReal<T>::type;
    const double tolerance_;
    int restart_;
    LinearWorkspace<T, Device> work_;
    LinearAlgebra<T, Device> algebra_;
    ct::Tensor krylov_;
    std::vector<std::vector<Wide>> h_;
    std::vector<std::vector<Wide>> g_;
    std::vector<std::vector<Wide>> sine_;
    std::vector<std::vector<double>> cosine_;
    std::vector<int> order_;
    std::vector<int> swaps_;
    int ld_ = 0;
    int dim_ = 0;
    int bands_ = 0;
    int stride_ = 0;
    int active_ = 0;

    T* vector(int slot);
    T* residual();
    T* solution();
    T* basis(int index);
    T* direction(int index);
    T* image(int index);
    bool start_cycle(const std::vector<double>& threshold, std::vector<T>* coefficients);
    void orthogonalize(int j, T* next, std::vector<T>* coefficients);
    bool update_qr(int b, int j, double norm);
    bool back_substitute(int b, int j, std::vector<Wide>* weights);
    void update_solution(int j,
                         const std::vector<int>& done,
                         const std::vector<std::vector<Wide>>& weights,
                         bool reconstruct,
                         std::vector<T>* coefficients);
    void retire(int band);
    void apply_swaps(int last);
    bool cycle(const LinearOperator<T, Device>& op,
               const LinearOperator<T, Device>& preconditioner,
               const std::vector<double>& threshold,
               bool reconstruct,
               int max_iterations,
               LinearSolveResult* result);

  public:
    LinearGMRES(double tolerance, const LinearSolveOptions& options, const diag_comm_info& comm);
    /** @brief Use the per-call budget and reconstruction policy; force_check requests an independent final application. */
    LinearSolveResult solve(const LinearOperator<T, Device>& op,
                            const LinearOperator<T, Device>& preconditioner,
                            int ld,
                            int nvec,
                            int dim,
                            T* x,
                            const T* b,
                            const T* initial_residual,
                            bool force_check,
                            const LinearSolveControl& control);
};
} // namespace hsolver
#endif
