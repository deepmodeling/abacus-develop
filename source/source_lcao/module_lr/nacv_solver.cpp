/** @file nacv_solver.cpp
 * @brief Certify CG solutions and reuse direct factorizations for the NAC response.
 */
#include "nacv_solver.h"
#include "nacv_tools.h"
#include "source_base/module_container/base/third_party/lapack.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/opt_cg.h"
#include "source_base/parallel_blkcyc.h"
#include "source_base/parallel_cell.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "utils/lr_util.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

// The diagnostic printer requires the stream/container/layout definitions above.
#include "utils/lr_util_print.h"

namespace LR
{
namespace
{
/** @brief Time the complete CG solve, including failed operator applications. */
class ZVectorCgTimer
{
  public:
    ZVectorCgTimer() { ModuleBase::timer::start("Z_vector", "solve_nac_cg"); }
    ~ZVectorCgTimer() { ModuleBase::timer::end("Z_vector", "solve_nac_cg"); }
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

void solve_nac_cg(
    double *const Z,
    double *R,
    const int &ld,
    const int &nstates,
    std::function<void(const double *const, double *const)> f_LZ,
    const bool test_force,
    const double tolerance,
    const int max_iterations)
{
    ModuleBase::TITLE("Z_vector", "solve_nac_cg");
    const ZVectorCgTimer timer;
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
    bool true_residual_ready = false;
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
                true_residual_ready = true;
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
    // A successful break already certified this unchanged Z. The update-limit
    // path still needs its own explicit operator application.
    if (!true_residual_ready)
    {
        f_LZ(Z, LP.data<double>());
        residual = response_residual(
            LP.data<double>(),
            R,
            correction.data(),
            size);
    }
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

void solve_nac_cg(
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

} // namespace LR

namespace LR
{
namespace
{
/** @brief Report numerical factorization/solve failures without accepting a response. */
void require_success(int info, const char *operation)
{
    if (info != 0)
    {
        throw std::runtime_error(std::string("NAC response ") + operation + " failed, info=" + std::to_string(info));
    }
}
} // namespace
NacvFactorization::NacvFactorization(std::vector<double> matrix, int dimension) : dimension_(dimension)
{
    const NacvTimerScope timer("NacvFactorization", "NacvFactorization");
    factors_ = std::move(matrix);
    pivots_.resize(dimension_);
    int info = 0;
    container::lapackConnector::getrf(
        dimension_,
        dimension_,
        factors_.data(),
        dimension_,
        pivots_.data(),
        info);
    require_success(info, "dgetrf");
}
void NacvFactorization::solve(std::vector<double> &rhs, int count)
{
    const NacvTimerScope timer("NacvFactorization", "solve");
    int info = 0;
    container::lapackConnector::getrs(
        'N',
        dimension_,
        count,
        factors_.data(),
        dimension_,
        pivots_.data(),
        rhs.data(),
        dimension_,
        info);
    require_success(info, "dgetrs");
}
#ifdef __MPI
NacvFactorization::NacvFactorization(std::vector<double> matrix, const Parallel_Blkcyclic &layout, bool cholesky)
    : layout_(&layout), cholesky_(cholesky)
{
    const NacvTimerScope timer("NacvFactorization", "NacvFactorization");
    dimension_ = layout.nrglb();
    factors_ = std::move(matrix);
    const std::size_t local_storage = std::max<std::size_t>(1, layout.size());
    factors_.resize(local_storage);
    if (cholesky_)
    {
        const std::size_t storage = std::max<std::size_t>(1, layout.size());
        std::vector<double> transposed(storage, 0.0);
        ScalapackConnector::tranc(
            dimension_,
            dimension_,
            1.0,
            factors_.data(),
            1,
            1,
            layout.desc(),
            0.0,
            transposed.data(),
            1,
            1,
            layout.desc());
        for (std::size_t i = 0; i < layout.size(); ++i)
        {
            factors_[i] = 0.5 * (factors_[i] + transposed[i]);
        }
        int descriptor[9];
        std::copy(layout.desc(), layout.desc() + 9, descriptor);
        const int info = ScalapackConnector::potrf(
            'U',
            dimension_,
            factors_.data(),
            descriptor);
        require_success(info, "pdpotrf (positive-definite Hessian required)");
    }
    else
    {
        const int pivot_count = layout.nrloc() + layout.nb();
        pivots_.resize(pivot_count);
        int info = 0;
        ScalapackConnector::getrf(
            dimension_,
            dimension_,
            factors_.data(),
            1,
            1,
            layout.desc(),
            pivots_.data(),
            &info);
        require_success(info, "pdgetrf");
    }
}
void NacvFactorization::solve(std::vector<double> &rhs, const Parallel_Blkcyclic &layout)
{
    const NacvTimerScope timer("NacvFactorization", "solve");
    const int count = layout.ncglb();
    int info = 0;
    if (cholesky_)
    {
        ScalapackConnector::potrs(
            'U',
            dimension_,
            count,
            factors_.data(),
            1,
            1,
            layout_->desc(),
            rhs.data(),
            1,
            1,
            layout.desc(),
            &info);
    }
    else
    {
        ScalapackConnector::getrs(
            'N',
            dimension_,
            count,
            factors_.data(),
            1,
            1,
            layout_->desc(),
            pivots_.data(),
            rhs.data(),
            1,
            1,
            layout.desc(),
            &info);
    }
    require_success(info, "distributed triangular solve");
}
#endif
} // namespace LR
