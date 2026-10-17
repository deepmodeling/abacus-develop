#ifndef ABACUS_LR_HSOLVER_LRTD_H
#define ABACUS_LR_HSOLVER_LRTD_H
#include "source_base/parallel_2d.h"
#include "source_base/macros.h"
#include <complex>
#include <string>
#include <vector>
namespace LR
{
template<typename T> class HamiltLR;
template<typename T> class HamiltULR;
template<typename T> using Real = typename GetTypeReal<T>::type;
namespace HSolver
{
/// Solve the distributed LR eigenproblem using the selected dense or iterative backend.
        template<typename T, typename THamilt>
        void solve(const THamilt& hm,
            T* psi,
            const int& dim, ///< local leading dimension (or nbasis)
            const int& nband,   ///< nstates in LR-TDDFT, not (nocc+nvirt)
            const int& nk,
            const std::vector<int>& nocc,
            const std::vector<int>& nvirt,
            const std::vector<Parallel_2D>& pX,
            double* eig,
            const std::string method,
            const Real<T>& diag_ethr, ///< threshold for diagonalization
            const std::vector<Real<T>>& precondition,
            const bool hermitian);
}
}
#endif
