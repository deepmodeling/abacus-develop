#include "zeqlin_solv.h"

#ifdef __MPI
#include <algorithm>
#include <cassert>
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

#include "source_base/module_external/scalapack_connector.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"

namespace LR
{
    /// A <- (A + A^H) / 2 on the 2D layout `pA`. The Cholesky solvers read only the upper
    /// triangle, so any numerical asymmetry of A would otherwise be dropped silently instead of
    /// averaged.
    template <typename T>
    void hermitize_2d(T* A, const Parallel_2D& pA)
    {
        const int n = pA.get_global_row_size();
        std::vector<T> AH(pA.get_local_size(), T(0));
        const T one(1.0);
        const T zero(0.0);
        ScalapackConnector::tranc(n, n, one, A, 1, 1, pA.desc, zero, AH.data(), 1, 1, pA.desc);
        for (std::size_t i = 0; i < AH.size(); ++i)
        {
            A[i] = 0.5 * (A[i] + AH[i]);
        }
    }

    template <typename T>
    void scalapack_linear_solver(T* A, T* B, const Parallel_2D& pA, const Parallel_2D& pB)
    {
        ModuleBase::TITLE("LR", "scalapack_linear_solver");
        ModuleBase::timer::start("LR", "scalapack_linear_solver");
        const int n = pA.get_global_row_size();
        const int nrhs = pB.get_global_col_size();
        assert(pA.get_global_col_size() == n);
        assert(pB.get_global_row_size() == n);
        assert(pA.blacs_ctxt == pB.blacs_ctxt);
        assert(pA.get_block_size() == pB.get_block_size());
        // p?gesv needs LOCr(M_A) + MB_A pivot entries
        std::vector<int> ipiv(pA.get_row_size() + pA.get_block_size());
        int info = 0;
        ScalapackConnector::gesv(n, nrhs, A, 1, 1, pA.desc, ipiv.data(), B, 1, 1, pB.desc, &info);
        if (info != 0)
        {
            ModuleBase::timer::end("LR", "scalapack_linear_solver");
            throw std::runtime_error("scalapack_linear_solver: p?gesv failed, info=" + std::to_string(info));
        }
        ModuleBase::timer::end("LR", "scalapack_linear_solver");
    }

    template <typename T>
    void scalapack_cholesky_linear_solver(T* A, T* B, const Parallel_2D& pA, const Parallel_2D& pB)
    {
        ModuleBase::TITLE("LR", "scalapack_cholesky_linear_solver");
        ModuleBase::timer::start("LR", "scalapack_cholesky_linear_solver");
        const int n = pA.get_global_row_size();
        const int nrhs = pB.get_global_col_size();
        assert(pA.get_global_col_size() == n);
        assert(pB.get_global_row_size() == n);
        assert(pA.blacs_ctxt == pB.blacs_ctxt);
        assert(pA.get_block_size() == pB.get_block_size());

        hermitize_2d(A, pA);

        // A = U^H U, U in the upper triangle. p?potrf's INFO is global output, so a
        // non-positive-definite A is reported identically on every rank.
        int desc_a[9];
        std::copy(pA.desc, pA.desc + 9, desc_a);
        int info = ScalapackConnector::potrf('U', n, A, desc_a);
        if (info != 0)
        {
            std::string message = "scalapack_cholesky_linear_solver: p?potrf failed, info="
                + std::to_string(info);
            if (info > 0)
            {
                message += " -- the matrix is not positive definite; "
                    "use lr_grad_solver=scalapack (LU) instead";
            }
            ModuleBase::timer::end("LR", "scalapack_cholesky_linear_solver");
            throw std::runtime_error(message);
        }
        ScalapackConnector::potrs('U', n, nrhs, A, 1, 1, pA.desc, B, 1, 1, pB.desc, &info);
        if (info != 0)
        {
            ModuleBase::timer::end("LR", "scalapack_cholesky_linear_solver");
            throw std::runtime_error("scalapack_cholesky_linear_solver: p?potrs failed, info="
                + std::to_string(info));
        }
        ModuleBase::timer::end("LR", "scalapack_cholesky_linear_solver");
    }

    template void scalapack_linear_solver<double>(double*, double*, const Parallel_2D&, const Parallel_2D&);
    template void scalapack_linear_solver<std::complex<double>>(std::complex<double>*, std::complex<double>*,
        const Parallel_2D&, const Parallel_2D&);
    template void scalapack_cholesky_linear_solver<double>(double*, double*, const Parallel_2D&,
        const Parallel_2D&);
    template void scalapack_cholesky_linear_solver<std::complex<double>>(std::complex<double>*,
        std::complex<double>*, const Parallel_2D&, const Parallel_2D&);

}
#endif
