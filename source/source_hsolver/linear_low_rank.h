#ifndef HSOLVER_LINEAR_LOW_RANK_H
#define HSOLVER_LINEAR_LOW_RANK_H
#include "source_hsolver/linear_algebra.h"
#include "source_hsolver/linear_operator.h"

namespace hsolver
{
/** @brief A bounded, paired-orthogonalized response history for one k point. */
template <typename T, typename Device>
class LinearResponse
{
  private:
    ct::Tensor directions_, images_;
    int rank_ = 0;

  public:
    int rank() const
    {
        return rank_;
    }
    void clear()
    {
        rank_ = 0;
    }
    const T* directions() const
    {
        return directions_.template data<T>();
    }
    const T* images() const
    {
        return images_.template data<T>();
    }

    /** @brief Replace the history only after a successful solve; discard unresolved directions. */
    void update(LinearAlgebra<T, Device>& algebra,
                int ld,
                int dim,
                int bands,
                const T* solution,
                const T* seed,
                const T* seed_residual,
                double tolerance,
                ct::Tensor* workspace);
};

/** @brief Diagonal inverse with an optional response or Galerkin coarse correction. */
template <typename T, typename Device>
class LinearLowRank final : public LinearOperator<T, Device>
{
  private:
    LinearAlgebra<T, Device>& algebra_;
    const T* diagonal_;
    const int dim_;
    const T* test_ = nullptr;
    int rank_ = 0;
    const LinearSmallLU* factor_ = nullptr;
    ct::Tensor correction_;

  public:
    LinearLowRank(LinearAlgebra<T, Device>& algebra, const T* diagonal, int dim) : algebra_(algebra), diagonal_(diagonal), dim_(dim)
    {
    }
    bool is_identity() const override
    {
        return diagonal_ == nullptr && rank_ == 0;
    }
    int rank() const
    {
        return rank_;
    }

    /** @brief Prepare a response correction with orthonormal images.
     *  @param images Borrowed projection vectors; keep them valid and unchanged until the last apply.
     */
    void prepare_response(int ld, int rank, const T* directions, const T* images);
    /** @brief Prepare a Galerkin coarse correction.
     *  @param basis Borrowed projection vectors; keep them valid and unchanged until the last apply.
     *  @param factor Borrowed factorization of basis^H*images; keep it valid and unchanged until the last apply.
     */
    void prepare_subspace(int ld, int rank, const T* basis, const T* images, const LinearSmallLU& factor);
    /** @brief Apply the fixed correction; throw LinearPreconditionerError if the coarse solve fails. */
    void apply(const T* x, T* y, int ld, int nvec) const override;

  private:
    /** @brief Form Z-D*W and retain the projection data for subsequent applications. */
    void prepare(int ld, int rank, const T* z, const T* w, const T* test, const LinearSmallLU* factor);
};
} // namespace hsolver
#endif
