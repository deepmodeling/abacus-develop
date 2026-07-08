#include <ATen/kernels/lapack.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace hsolver {
namespace {

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
        std::vector<Scalar> r(b, b + n * n);
        potrf(n, r.data());
        trtri(n, r.data());

        std::vector<Scalar> c(n * n, Scalar(0));
        for (int j = 0; j < n; ++j)
        {
            for (int i = 0; i < n; ++i)
            {
                Scalar sum = Scalar(0);
                for (int p = 0; p <= i; ++p)
                {
                    const Scalar rip = r[p + i * n];
                    if (rip == Scalar(0))
                        continue;
                    for (int q = 0; q <= j; ++q)
                    {
                        const Scalar rqj = r[q + j * n];
                        if (rqj != Scalar(0))
                            sum += std::conj(rip) * a[p + q * n] * rqj;
                    }
                }
                c[i + j * n] = sum;
            }
        }
        for (const Scalar& cij : c)
        {
            if (!std::isfinite(std::real(cij))
                || !std::isfinite(std::imag(cij)))
                throw std::runtime_error("PPCG: reduced matrix is non-finite.");
        }

        syevd(n, c.data(), w);
        for (int j = 0; j < n; ++j)
        {
            if (!std::isfinite(w[j]))
                throw std::runtime_error("PPCG: syevd returned non-finite eigenvalue.");
        }

        std::fill(a, a + n * n, Scalar(0));
        for (int j = 0; j < n; ++j)
        {
            Real nrm2 = 0;
            for (int i = 0; i < n; ++i)
            {
                Scalar sum = Scalar(0);
                for (int p = i; p < n; ++p)
                    sum += r[i + p * n] * c[p + j * n];
                if (!std::isfinite(std::real(sum))
                    || !std::isfinite(std::imag(sum)))
                    throw std::runtime_error("PPCG: back-transformed eigenvector is non-finite.");
                a[i + j * n] = sum;
                nrm2 += static_cast<Real>(std::norm(sum));
            }
            if (nrm2 <= static_cast<Real>(1e-30))
                throw std::runtime_error("PPCG: back-transformed eigenvector is zero.");
        }
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

} // anonymous namespace
} // namespace hsolver
