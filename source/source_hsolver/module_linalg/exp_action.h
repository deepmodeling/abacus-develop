#ifndef HSOLVER_EXP_ACTION_H
#define HSOLVER_EXP_ACTION_H

#include "source_hsolver/hs_operator.h"
#include "source_hsolver/module_linalg/linear_algebra.h"
#include "source_hsolver/module_linalg/orb_buffer.h"

#include <memory>

namespace hsolver
{
enum class ExpMethod
{
    taylor,
    lanczos
};

/** @brief Explicit controls for a frozen Hermitian Hamiltonian exponential. */
struct ExpControls
{
    ExpMethod expansion = ExpMethod::lanczos;
    int taylor_order = 4;
    double rel_tol = 0.0;
    int krylov_max = 32;
};

/** @brief Outcome of an exponential evaluation; retry requests a shorter interval. */
enum class ExpState
{
    success,
    retry,
    failed
};

/** @brief Error is relative to the input norm and only estimated by Lanczos. */
struct ExpResult
{
    ExpState state = ExpState::success;
    bool has_err_est = false;
    double rel_err_est = 0.0;
    int order_or_dim = 0;
    int bisect_count = 0;
    int reorth_count = 0;
    int failed_band = -1;
    long long hpsi_calls = 0;
    long long hpsi_cols = 0;
};

/** @brief Apply exp(-i*alpha*H_Ry); caller selects k, buffers must not alias. */
template <typename T, typename Device>
class ExpAction
{
  protected:
    ExpControls exp_ctrl_;
    LinearAlgebra<T, Device> algebra_;

  public:
    ExpAction(const ExpControls& exp_ctrl, const diag_comm_info& comm);
    virtual ~ExpAction() = default;
    /** @brief Effective relative tolerance after resolving the precision-dependent default. */
    double rel_tol() const
    {
        return exp_ctrl_.rel_tol;
    }
    virtual ExpResult apply(const HSOperator<T, Device>& op, const OrbLayout& layout, double alpha, const T* input, T* output) = 0;
    virtual ExpResult apply_pair(const HSOperator<T, Device>& op, const OrbLayout& layout, double alpha, const T* input, T* half, T* full)
        = 0;
};

/** @brief Construct a backend without reading global input parameters. */
template <typename T, typename Device>
std::unique_ptr<ExpAction<T, Device>> make_exp_action(const ExpControls& exp_ctrl, const diag_comm_info& comm);
} // namespace hsolver
#endif
