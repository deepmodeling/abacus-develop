#include "nacv_factorization.h"
#include "nacv_timer.h"
#include "source_base/module_container/base/third_party/lapack.h"
#ifdef __MPI
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/parallel_blkcyc.h"
#endif
#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>
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
