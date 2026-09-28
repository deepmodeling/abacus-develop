#ifndef HSOLVER_LINEAR_SOLVER_TYPES_H
#define HSOLVER_LINEAR_SOLVER_TYPES_H

#include <cstdint>

namespace hsolver
{

enum class LinearSolveStatus
{
    converged,
    breakdown,
    max_iterations,
    residual_mismatch
};

/** @brief Result checked against b - A*x in the original column order. */
struct LinearSolveResult
{
    LinearSolveStatus status = LinearSolveStatus::max_iterations;
    double max_residual = 0.0;
    int iterations = 0;
    int failed_band = -1;
    int restarts = 0;
    std::int64_t operator_calls = 0;
    std::int64_t operator_columns = 0;
};

enum class LinearMethod
{
    bicgstab,
    cgs
};

/** @brief Numerical policy; tolerance zero selects the precision-aware default. */
struct LinearSolveOptions
{
    LinearMethod method = LinearMethod::bicgstab;
    double tolerance = 0.0;
    int max_iterations = 500;
};

const char* linear_status_name(const LinearSolveStatus status);

} // namespace hsolver
#endif
