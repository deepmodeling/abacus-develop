#ifndef HSOLVER_CN_SUBSPACE_H
#define HSOLVER_CN_SUBSPACE_H
#include "source_hsolver/linear_algebra.h"

namespace hsolver
{
/** @brief CN projection derived from A*U = 2*U - B without another Hamiltonian application. */
template <typename T, typename Device>
class CNSubspace
{
  private:
    ct::Tensor image_;
    ct::Tensor seed_;
    ct::Tensor residual_;
    LinearSmallLU factor_;

  public:
    bool prepare(LinearAlgebra<T, Device>& algebra, int ld, int dim, int bands, const T* u, const T* b);
    T* image()
    {
        return image_.template data<T>();
    }
    T* seed()
    {
        return seed_.template data<T>();
    }
    T* residual()
    {
        return residual_.template data<T>();
    }
    const LinearSmallLU& factor() const
    {
        return factor_;
    }
};
} // namespace hsolver
#endif
