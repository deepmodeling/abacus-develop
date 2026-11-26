#ifndef ABACUS_LR_NACV_SOLVER_H
#define ABACUS_LR_NACV_SOLVER_H

#include <complex>
#include <functional>
#include <vector>

namespace LR
{
/**
 * @brief Solve a distributed real orbital-response equation and certify its residual.
 * @param[out] Z Local state-major solution, with nstates blocks of ld elements.
 * @param[in] R Local right-hand side, in Rydberg for the physical A+B response operator.
 * @param[in] ld Local vector length, including zero on ranks without entries.
 * @param[in] nstates Number of right-hand sides solved with a combined residual norm.
 * @param[in] f_LZ Collective linear operator; writes all local output entries.
 * @param[in] test_force Whether to print the full final vector for diagnostics.
 * @param[in] tolerance Positive finite absolute global Euclidean norm bound on LZ-R.
 * The physical response residual and tolerance have units Ry; Z is dimensionless.
 * @param[in] max_iterations Positive total update limit, including residual corrections.
 * @pre Dimensions are the prepared NAC response layout; controls passed INPUT validation.
 * @throw std::runtime_error If the explicitly recomputed residual is not below tolerance.
 * @details A converged recurrence is checked against LZ-R. When these differ,
 * retain Z and restart a correction equation with the true residual, using the
 * same tolerance and remaining global update budget. Apply a correction update
 * before another recursive-norm check, so rounding differences cannot cause
 * an unlimited sequence of zero-update restarts. R is never modified.
 */
void solve_nac_cg(
    double *const Z,
    double *R,
    const int &ld,
    const int &nstates,
    std::function<void(const double *const, double *const)> f_LZ,
    const bool test_force,
    const double tolerance,
    const int max_iterations);

/**
 * @brief Reject the unavailable complex CG response solver explicitly.
 * @param[out] Z Reserved local solution storage; not modified.
 * @param[in] R Reserved local right-hand side; not modified.
 * @param[in] ld Local vector length.
 * @param[in] nstates Number of right-hand sides.
 * @param[in] f_LZ Reserved complex operator.
 * @param[in] test_force Reserved diagnostic switch.
 * @param[in] tolerance Reserved absolute residual tolerance.
 * @param[in] max_iterations Reserved iteration bound.
 * @throw std::runtime_error Always: complex CG response is not implemented.
 */
void solve_nac_cg(
    std::complex<double> *const Z,
    std::complex<double> *R,
    const int &ld,
    const int &nstates,
    std::function<void(const std::complex<double> *const, std::complex<double> *const)> f_LZ,
    const bool test_force,
    const double tolerance,
    const int max_iterations);

} // namespace LR

class Parallel_Blkcyclic;
namespace LR
{
/**
 * @brief Own one real NAC response factorization for repeated right-hand sides.
 * @details Matrices are column-major. Distributed layouts are borrowed and
 * must outlive the factorization. No matrix inverse or tolerance change is used.
 */
class NacvFactorization
{
  public:
    /**
     * @param[in] matrix Complete orbital Hessian in Ry; ownership is transferred.
     * @param[in] dimension Global matrix dimension.
     */
    NacvFactorization(std::vector<double> matrix, int dimension);
#ifdef __MPI
    /**
     * @param[in] matrix Local Hessian block in Ry; ownership is transferred.
     * @param[in] layout Square block-cyclic layout, with the existing BLACS context.
     * @param[in] cholesky Use the existing symmetrized positive-definite algorithm.
     */
    NacvFactorization(std::vector<double> matrix, const Parallel_Blkcyclic &layout, bool cholesky);
    /**
     * @param[in,out] rhs Local RHS blocks in Ry, overwritten with dimensionless Z.
     * @param[in] layout RHS layout with the same context and block size as the Hessian.
     */
    void solve(std::vector<double> &rhs, const Parallel_Blkcyclic &layout);
#endif
    /**
     * @param[in,out] rhs Complete RHS columns in Ry, overwritten with dimensionless Z.
     * @param[in] count Number of RHS columns.
     */
    void solve(std::vector<double> &rhs, int count);

  private:
    int dimension_ = 0;
    const Parallel_Blkcyclic *layout_ = nullptr;
    bool cholesky_ = false;
    std::vector<double> factors_;
    std::vector<int> pivots_;
};
} // namespace LR

#endif
