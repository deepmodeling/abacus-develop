#ifndef HSOLVER_LINEAR_WORKSPACE_H
#define HSOLVER_LINEAR_WORKSPACE_H

#include "source_base/macros.h"
#include "source_base/module_container/ATen/core/tensor.h"
#include "source_hsolver/diag_comm_info.h"
#include "source_hsolver/linear_operator.h"
#include "source_hsolver/linear_solver_types.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace hsolver
{

/** @brief Balance solve timing on both normal and exceptional exits. */
class LinearSolveTimer
{
  private:
    const char* solver_;

  public:
    explicit LinearSolveTimer(const char* solver);
    ~LinearSolveTimer();
    LinearSolveTimer(const LinearSolveTimer&) = delete;
    LinearSolveTimer& operator=(const LinearSolveTimer&) = delete;
};

/** @brief Reusable device arrays and collective block algebra for linear solvers. */
template <typename T, typename Device>
class LinearWorkspace
{
  private:
    ct::Tensor vectors_;
    ct::Tensor dots_;
    ct::Tensor partial_;
    ct::Tensor coefficients_;
    ct::Tensor mask_;
    ct::Tensor permutation_;
    std::vector<T> host_dots_;
    std::vector<T> host_coefficients_;
    int capacity_ = 0;
    std::int64_t operator_calls_ = 0;
    std::int64_t operator_columns_ = 0;
    const diag_comm_info comm_;

  public:
    explicit LinearWorkspace(const diag_comm_info& comm);
    void prepare(const int ld, const int dim, const int nvec, const T* x, const T* b, const double tolerance, const int max_iter);
    /** @brief Initialize recurrence storage once per cycle, including padding. */
    void clear();
    void reset_statistics();
    void apply(const LinearOperator<T, Device>& op, const T* x, T* y, int ld, int nvec);
    void cgs_direction(int ld, int dim, int nvec, T* p, T* u, const T* r, const T* q, const T* beta);
    void cgs_alpha(int ld, int dim, int nvec, T* q, T* direction, const T* u, const T* v, const T* alpha);
    /** @brief Finish the CGS step using the coefficients uploaded by cgs_alpha. */
    void cgs_finish(int ld, int dim, int nvec, T* x, T* r, const T* direction, const T* ad);
    T* data(const int slot);
    void dot(const int ld, const int dim, const int nvec, const T* x, const T* y, T* out);
    /** @brief Return two packed global inner products in a single collective. */
    void dot_pair(int ld,int dim,int nvec,const T* x,const T* y,const T* z,const T* w,T* first,T* second);
    /** @brief Upload per-column coefficients and update all unmasked columns. */
    void batch(int ld,int dim,int nvec,T* out,const T* x,const T* y,T a,T b,
               const T* ca,const T* cb,const int* skip);
    /** @brief Reorder solver columns after host convergence decisions. */
    void swap_columns(int ld,int dim,const std::vector<int>& pairs);
    /** @brief Restore original column order into caller storage. */
    void restore(int ld,int dim,int nvec,const std::vector<int>& order,const T* source,T* destination);
    /** @brief Upload two coefficient arrays and execute a fused recurrence. */
    void bicg_update(int ld,int dim,int nvec,bool direction,T* out,const T* x,const T* y,const T* a,const T* b);
    void copy(const int ld, const int dim, const int nvec, const T* x, T* y) const;
    void residual(const LinearOperator<T, Device>& op, const int ld, const int dim, const int nvec, const T* x, const T* b, T* r);
    void verify(const LinearOperator<T, Device>& op,
                const int ld,
                const int dim,
                const int nvec,
                const T* x,
                const T* b,
                const std::vector<double>& threshold,
                T* scratch,
                LinearSolveResult* result);
};

template <typename T>
bool linear_bad_divisor(const T value)
{
    using Real = typename GetTypeReal<T>::type;
    return !std::isfinite(std::abs(value)) || std::abs(value) <= std::numeric_limits<Real>::min();
}

template <typename T>
double linear_norm(const T value)
{
    return std::isfinite(std::abs(value)) ? std::sqrt(std::max(0.0, static_cast<double>(std::real(value))))
                                          : std::numeric_limits<double>::infinity();
}

} // namespace hsolver
#endif
