#include "zeqlin_solv.h"

#ifdef __MPI
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

#include "source_base/module_external/blacs_connector.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef __ELPA
#include "source_hsolver/module_genelpa/elpa_new.h"
#ifdef I // avoid conflict with the macro defined by ELPA
#undef I
#endif
#endif

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
        for (std::size_t i = 0; i < AH.size(); ++i) { A[i] = 0.5 * (A[i] + AH[i]); }
    }

#ifdef __ELPA
    /// Number of diagonal entries of the distributed Cholesky factor U that are not real, positive
    /// and finite, summed over the BLACS grid of `pA`. A failed potrf leaves its non-positive pivot
    /// on the diagonal, so this is nonzero exactly when the factorization broke down.
    template <typename T>
    int count_bad_cholesky_pivots(const T* U, const Parallel_2D& pA)
    {
        int nbad = 0;
        for (int lc = 0; lc < pA.get_col_size(); ++lc)
        {
            const int lr = pA.global2local_row(pA.local2global_col(lc));
            if (lr < 0) { continue; }
            const T u = U[static_cast<std::size_t>(lc) * pA.get_row_size() + lr];
            const double re = std::real(u);
            const bool good = std::isfinite(re) && re > 0.0 && std::abs(std::imag(u)) <= 1e-12 * re;
            if (!good) { ++nbad; }
        }
        char scope[] = "All";
        char top[] = " ";
        Cigsum2d(pA.blacs_ctxt, scope, top, 1, 1, &nbad, 1, -1, -1);
        return nbad;
    }
#endif

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

        // A = U^H U, U in the upper triangle. Unlike ELPA's, p?potrf's INFO is global output, so a
        // non-positive-definite A is reported identically on every rank.
        int desc_a[9];
        std::copy(pA.desc, pA.desc + 9, desc_a);
        int info = ScalapackConnector::potrf('U', n, A, desc_a);
        if (info != 0)
        {
            throw std::runtime_error("scalapack_cholesky_linear_solver: p?potrf failed, info="
                + std::to_string(info) + " -- the matrix is not positive definite; "
                "use the LU-based 'scalapack' solver instead");
        }
        ScalapackConnector::potrs('U', n, nrhs, A, 1, 1, pA.desc, B, 1, 1, pB.desc, &info);
        if (info != 0)
        {
            throw std::runtime_error("scalapack_cholesky_linear_solver: p?potrs failed, info="
                + std::to_string(info));
        }
        ModuleBase::timer::end("LR", "scalapack_cholesky_linear_solver");
    }

    template <typename T>
    void elpa_linear_solver(T* A, T* B, const Parallel_2D& pA, const Parallel_2D& pB)
    {
        ModuleBase::TITLE("LR", "elpa_linear_solver");
        ModuleBase::timer::start("LR", "elpa_linear_solver");
#ifdef __ELPA
        const int n = pA.get_global_row_size();
        const int nrhs = pB.get_global_col_size();
        assert(pA.get_global_col_size() == n);
        assert(pB.get_global_row_size() == n);
        assert(pA.blacs_ctxt == pB.blacs_ctxt);
        assert(pA.get_block_size() == pB.get_block_size());

        hermitize_2d(A, pA);

        int status = 0;
        if (elpa_init(20210430) != ELPA_OK)
        {
            throw std::runtime_error("elpa_linear_solver: ELPA API version not supported");
        }
        elpa_t handle = elpa_allocate(&status);
        if (status != ELPA_OK) { throw std::runtime_error("elpa_linear_solver: elpa_allocate failed"); }
        elpa_set(handle, "na", n, &status);
        elpa_set(handle, "nev", n, &status);
        elpa_set(handle, "local_nrows", pA.get_row_size(), &status);
        elpa_set(handle, "local_ncols", pA.get_col_size(), &status);
        elpa_set(handle, "nblk", pA.get_block_size(), &status);
        elpa_set(handle, "mpi_comm_parent", MPI_Comm_c2f(pA.comm()), &status);
        elpa_set(handle, "process_row", pA.get_coord_row(), &status);
        elpa_set(handle, "process_col", pA.get_coord_col(), &status);
#ifdef _OPENMP
        const int num_threads = omp_get_max_threads();
#else
        const int num_threads = 1;
#endif
        elpa_set(handle, "omp_threads", num_threads, &status);
        if (elpa_setup(handle) != ELPA_OK) { throw std::runtime_error("elpa_linear_solver: elpa_setup failed"); }

        // A = U^H U, U in the upper triangle.
        // The returned status cannot be trusted: ELPA 2022.11 declares the C binding's `error` as
        // intent(in), so a failed Cholesky still reports ELPA_OK. Check the pivots instead.
        // On a non-positive-definite A only the rank owning the failing diagonal block returns
        // (elpa_cholesky_template.F90) while the others wait in a broadcast, so with several
        // ranks the check below is never reached and the run hangs.
        elpa_cholesky(handle, A, &status);
        // No `elpa_uninit` here: the KS solver may still hold ELPA handles of its own.
        elpa_deallocate(handle, &status);
        if (count_bad_cholesky_pivots(A, pA) > 0)
        {
            throw std::runtime_error("elpa_linear_solver: ELPA Cholesky failed -- the matrix is not "
                "positive definite; use the LU-based 'scalapack' solver instead");
        }

        int info = 0;
        ScalapackConnector::potrs('U', n, nrhs, A, 1, 1, pA.desc, B, 1, 1, pB.desc, &info);
        if (info != 0)
        {
            throw std::runtime_error("elpa_linear_solver: p?potrs failed, info=" + std::to_string(info));
        }
#else
        throw std::runtime_error("elpa_linear_solver: ABACUS was built without ELPA; rebuild with "
            "ENABLE_ELPA=ON or use the 'scalapack' solver");
#endif
        ModuleBase::timer::end("LR", "elpa_linear_solver");
    }

    template void scalapack_linear_solver<double>(double*, double*, const Parallel_2D&, const Parallel_2D&);
    template void scalapack_linear_solver<std::complex<double>>(std::complex<double>*, std::complex<double>*,
        const Parallel_2D&, const Parallel_2D&);
    template void scalapack_cholesky_linear_solver<double>(double*, double*, const Parallel_2D&,
        const Parallel_2D&);
    template void scalapack_cholesky_linear_solver<std::complex<double>>(std::complex<double>*,
        std::complex<double>*, const Parallel_2D&, const Parallel_2D&);
    template void elpa_linear_solver<double>(double*, double*, const Parallel_2D&, const Parallel_2D&);
    template void elpa_linear_solver<std::complex<double>>(std::complex<double>*, std::complex<double>*,
        const Parallel_2D&, const Parallel_2D&);
}
#endif
