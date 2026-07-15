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

    static void sygvd(int n, Scalar* a, Scalar* b, Real* w)
    {
        std::vector<Scalar> eigenvectors(n * n);
        container::kernels::lapack_hegvd<Scalar, Device>()(
            n, n, a, b, w, eigenvectors.data());
        std::copy(eigenvectors.begin(), eigenvectors.end(), a);
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
