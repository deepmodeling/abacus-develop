#ifndef ABACUS_LR_NACV_FACTORIZATION_H
#define ABACUS_LR_NACV_FACTORIZATION_H
#include <vector>
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
