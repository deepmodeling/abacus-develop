#include <ATen/kernels/lapack.h>

#include <cstdlib>
#include <fstream>

namespace hsolver {

// =============================================================================
// LAPACK wrapper (specialized per real type)
// =============================================================================
namespace {

template <typename T, typename Real>
Real max_generalized_residual(
    const T* hpsi,
    const T* spsi,
    const Real* eigenvalue,
    int ld,
    int n_dim,
    int ncol)
{
    Real max_res = 0;
    for (int j = 0; j < ncol; ++j)
    {
        Real nrm2 = 0;
        for (int ig = 0; ig < n_dim; ++ig)
        {
            const T r = hpsi[ig + j * ld] - T(eigenvalue[j]) * spsi[ig + j * ld];
            nrm2 += static_cast<Real>(std::norm(r));
        }
        max_res = std::max(max_res, std::sqrt(nrm2));
    }
    return max_res;
}

template <typename Scalar>
struct HermitianLapack
{
    using Real = typename container::GetTypeReal<Scalar>::type;
    using Device = container::DEVICE_CPU;

    static void syevd(int n, Scalar* a, Real* w)
    {
        container::kernels::lapack_heevd<Scalar, Device>()(n, a, n, w);
    }

    static void sygvd(int n, Scalar* a, Scalar* b, Real* w)
    {
        std::vector<Scalar> eigvec(n * n, Scalar(0));
        container::kernels::lapack_hegvd<Scalar, Device>()(n, n, a, b, w, eigvec.data());
        for (int j = 0; j < n; ++j)
        {
            if (!std::isfinite(w[j]))
                throw std::runtime_error("PPCG: hegvd returned non-finite eigenvalue.");

            Real nrm2 = 0;
            for (int i = 0; i < n; ++i)
                nrm2 += static_cast<Real>(std::norm(eigvec[i + j * n]));
            if (nrm2 <= static_cast<Real>(1e-30))
                throw std::runtime_error("PPCG: hegvd returned a zero eigenvector.");
        }
        std::copy(eigvec.begin(), eigvec.end(), a);
    }

    static void potrf(int n, Scalar* a)
    {
        Real diag_max = 0;
        for (int i = 0; i < n; ++i)
            diag_max = std::max(diag_max, std::abs(a[i + i * n]));
        std::vector<Scalar> a0(a, a + n * n);

        for (const Real shift : {Real(0), Real(1e-12), Real(1e-10), Real(1e-8),
                                 Real(1e-6), Real(1e-4), Real(1e-3), Real(1e-2),
                                 Real(1e-1), Real(1)})
        {
            std::copy(a0.begin(), a0.end(), a);
            if (shift > 0)
            {
                for (int i = 0; i < n; ++i)
                    a[i + i * n] += Scalar(shift * std::max(diag_max, Real(1)), 0);
            }
            try
            {
                container::kernels::lapack_potrf<Scalar, Device>()('U', n, a, n);
                return;
            }
            catch (const std::runtime_error&)
            {
                // Try the next diagonal shift.
            }
        }
        throw std::runtime_error("PPCG: potrf failed.");
    }

    static void trtri(int n, Scalar* a)
    {
        container::kernels::lapack_trtri<Scalar, Device>()('U', 'N', n, a, n);
    }
};

template <typename T>
inline void set_zero(std::vector<T>& x)
{
    std::fill(x.begin(), x.end(), T(0));
}

} // anonymous namespace

} // namespace hsolver
