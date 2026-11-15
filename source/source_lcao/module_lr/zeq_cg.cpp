#include "zeq_cg.h"
#include "source_base/opt_cg.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "utils/lr_util.h"
#include "utils/lr_util_print.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>
namespace LR
{
namespace
{
/** @brief Time the complete CG solve, including failed operator applications. */
class ZVectorCgTimer
{
  public:
    ZVectorCgTimer() { ModuleBase::timer::start("Z_vector", "solve_Z_CG"); }
    ~ZVectorCgTimer() { ModuleBase::timer::end("Z_vector", "solve_Z_CG"); }
};

/**
 * @brief Form the correction RHS from the explicitly reapplied linear operator.
 * @param[in] product Local LZ values for the accumulated solution.
 * @param[in] rhs Original local right-hand side, unchanged throughout the solve.
 * @param[out] correction Local R-LZ values, with the same distributed layout.
 * @param[in] size Number of local elements; zero-storage ranks still participate.
 * @return Global Euclidean norm of R-LZ in the response equation's Ry units.
 */
double response_residual(
    const double *product,
    const double *rhs,
    double *correction,
    int size)
{
    double squared = 0.0;
    for (int index = 0; index < size; ++index)
    {
        correction[index] = rhs[index] - product[index];
        squared += correction[index] * correction[index];
    }
    Parallel_Reduce::reduce_all(squared);
    return std::sqrt(squared);
}
} // namespace

void solve_Z_CG(
    double *const Z,
    double *R,
    const int &ld,
    const int &nstates,
    std::function<void(const double *const, double *const)> f_LZ,
    const bool test_force,
    const double tolerance,
    const int max_iterations)
{
    ModuleBase::TITLE("Z_vector", "solve_Z_CG");
    const ZVectorCgTimer timer;
    if (!std::isfinite(tolerance) || tolerance <= 0.0 || max_iterations < 1 || ld < 0 || nstates < 1 ||
        ld > std::numeric_limits<int>::max() / nstates)
    {
        throw std::invalid_argument(
            "Z-vector CG requires valid local dimensions, a positive finite tolerance and positive iteration bound");
    }
    const int maxiter = max_iterations;
    const double tol = tolerance;
    const int size = nstates * ld;
    container::Tensor P = LR_Util::newTensor<double>({nstates, ld});
    container::Tensor LP = LR_Util::newTensor<double>({nstates, ld});
    std::vector<double> correction(size);
    std::fill_n(Z, size, 0.0);
    ModuleBase::Opt_CG cg;
    cg.allocate(size);
    cg.init_b(R);
    std::cout << "Start solving Z-vector equaiton with CG method ..." << std::endl;
    int iterations = 0;
    int replacements = 0;
    bool correction_update_pending = false;
    double residual = std::numeric_limits<double>::infinity();
    while (iterations < maxiter)
    {
        cg.next_direct(LP.data<double>(), 0, P.data<double>());
        residual = cg.get_residual();
        std::cout << "iter=" << iterations << " residual=" << residual << std::endl;
        if (!std::isfinite(residual))
        {
            throw std::runtime_error("Z-vector CG produced a nonfinite recursive residual");
        }
        if (residual < tol && !correction_update_pending)
        {
            f_LZ(Z, LP.data<double>());
            residual = response_residual(
                LP.data<double>(),
                R,
                correction.data(),
                size);
            if (residual < tol)
            {
                break;
            }
            if (!std::isfinite(residual))
            {
                throw std::runtime_error("Z-vector CG produced a nonfinite true residual");
            }
            // Keep the accumulated Z and solve L(delta Z)=R-LZ. Refreshing the
            // recurrence replaces its drifted residual without accepting it or
            // discarding previous updates. The original global update limit holds.
            ++replacements;
            std::cout << "Z-vector CG residual replacement=" << replacements << " true residual=" << residual
                      << std::endl;
            cg.refresh(0, correction.data());
            std::fill_n(LP.data<double>(), size, 0.0);
            // The BLAS recurrence norm can round below tol while the explicit
            // norm equals it. Apply this correction direction before testing
            // the recurrence again, so a replacement cannot restart forever
            // without spending an update from the configured budget.
            correction_update_pending = true;
            continue;
        }
        f_LZ(P.data<double>(), LP.data<double>());
        int ifPD = 0;
        const double step = cg.step_length(LP.data<double>(), P.data<double>(), ifPD);
        const double *direction = P.data<double>();
        for (int i = 0; i < size; ++i)
        {
            Z[i] += step * direction[i];
        }
        ++iterations;
        correction_update_pending = false;
    }
    // Certification also covers termination at the configured update limit.
    f_LZ(Z, LP.data<double>());
    residual = response_residual(
        LP.data<double>(),
        R,
        correction.data(),
        size);
    std::cout << "Z-vector CG true residual=" << residual << " tolerance=" << tol << std::endl;
    if (!(residual < tol))
    {
        std::ostringstream message;
        message << "Z-vector CG did not converge after " << iterations << " iterations (true residual=" << residual
                << ", tolerance=" << tol << "). Excited-state forces are not valid.";
        throw std::runtime_error(message.str());
    }
    if (test_force)
    {
        std::cout << "Final Z-vector:" << std::endl;
        LR_Util::print_value(Z, nstates, ld);
    }
}

void solve_Z_CG(
    std::complex<double> *const Z,
    std::complex<double> *R,
    const int &ld,
    const int &nstates,
    std::function<void(const std::complex<double> *const, std::complex<double> *const)> f_LZ,
    const bool test_force,
    const double tolerance,
    const int max_iterations)
{
    throw std::runtime_error("complex Z-vector solver is not implemented yet");
}

void solve_Z_CG(
    double *Z,
    double *R,
    const int &ld,
    const int &nstates,
    std::function<void(const double *, double *)> f_LZ,
    bool test_force)
{
    const double tolerance = 1e-6;
    const int max_iterations = 100;
    solve_Z_CG(
        Z,
        R,
        ld,
        nstates,
        f_LZ,
        test_force,
        tolerance,
        max_iterations);
}

void solve_Z_CG(
    std::complex<double> *Z,
    std::complex<double> *R,
    const int &ld,
    const int &nstates,
    std::function<void(const std::complex<double> *, std::complex<double> *)> f_LZ,
    bool test_force)
{
    const double tolerance = 1e-6;
    const int max_iterations = 100;
    solve_Z_CG(
        Z,
        R,
        ld,
        nstates,
        f_LZ,
        test_force,
        tolerance,
        max_iterations);
}
} // namespace LR
