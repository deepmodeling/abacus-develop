#include "source_pw/module_pwdft/kernels/force_op.h"

#include <gtest/gtest.h>
#include <algorithm>
#include <complex>
#include <vector>
#ifdef __CUDA
#include <cuda_runtime.h>
#endif

namespace
{
template <typename Real>
struct DeltaSpinFixture
{
    int nh[2] = {2, 1};
    int na[2] = {1, 1};
    int isk[2] = {0, 1};
    const int nbands = 3;
    const int occupied = 2;
    const int nkb = 3;
    const int wg_nc = 4;
    const int ik = 1;
    Real tpiba = 1.7;
    std::vector<Real> weights = {0.9, 0.8, 0.7, 0.6, 0.3, 0.17, 0.0, 0.0};
    std::vector<Real> lambda;
    std::vector<std::complex<Real>> b;
    std::vector<std::complex<Real>> db;

    void prepare(const int npol, const int axis)
    {
        lambda.assign(6, 0);
        lambda[axis] = 0.7;
        lambda[3 + axis] = -0.4;
        b.resize(nbands * npol * nkb);
        db.resize(3 * b.size());
        for (std::size_t i = 0; i < b.size(); ++i)
        {
            b[i] = std::complex<Real>(0.13 + 0.017 * i, -0.21 + 0.009 * i);
        }
        for (std::size_t i = 0; i < db.size(); ++i)
        {
            db[i] = std::complex<Real>(-0.07 + 0.003 * i, 0.11 - 0.005 * i);
        }
    }

    // Independent oracle: form physical Pauli moments of b(R), then -d(lambda.M)/dR.
    // Promote the input to double before differencing, including for float kernels.
    double energy(const int atom, const int dir, const int npol, const double displacement) const
    {
        double result = 0;
        const int first = atom == 0 ? 0 : 2;
        const int count = nh[atom];
        const int sign = isk[ik] == 0 ? 1 : -1;
        for (int ib = 0; ib < occupied; ++ib)
        {
            for (int ip = first; ip < first + count; ++ip)
            {
                const int index = ib * npol * nkb + ip;
                const int derivative = dir * nbands * npol * nkb + index;
                const std::complex<double> up = std::complex<double>(b[index])
                    + displacement * double(tpiba) * std::complex<double>(db[derivative]);
                double value = sign * lambda[3 * atom + 2] * std::norm(up);
                if (npol == 2)
                {
                    const std::complex<double> down = std::complex<double>(b[index + nkb])
                        + displacement * double(tpiba) * std::complex<double>(db[derivative + nkb]);
                    const std::complex<double> cross = std::conj(up) * down;
                    value = lambda[3 * atom] * 2 * cross.real()
                        + lambda[3 * atom + 1] * 2 * cross.imag()
                        + lambda[3 * atom + 2] * (std::norm(up) - std::norm(down));
                }
                result += weights[ik * wg_nc + ib] * value;
            }
        }
        return result;
    }

    void check(const std::vector<Real>& force, const int npol, const double tolerance) const
    {
        const double h = 1e-4;
        for (int atom = 0; atom < 2; ++atom)
        {
            for (int dir = 0; dir < 3; ++dir)
            {
                const double ep = energy(atom, dir, npol, h);
                const double em = energy(atom, dir, npol, -h);
                const double expected = -(ep - em) / (2 * h);
                // Start with a nonzero force to also verify additive accumulation.
                EXPECT_NEAR(force[3 * atom + dir] - Real(0.125), expected, tolerance)
                    << "atom=" << atom << " dir=" << dir << " npol=" << npol;
            }
        }
    }

    void cpu(const int npol, std::vector<Real>& force) const
    {
        hamilt::cal_force_nl_op<Real, base_device::DEVICE_CPU>()(
            nullptr, occupied, wg_nc, 2, 3, nbands, ik, nkb, npol, nh, na, tpiba,
            weights.data(), lambda.data(), isk, b.data(), db.data(), force.data());
    }

#ifdef __CUDA
    void gpu(const int npol, std::vector<Real>& force) const
    {
        // Managed buffers keep the test independent of ABACUS memory helpers.
        int* metadata = nullptr;
        Real* scalars = nullptr;
        std::complex<Real>* overlaps = nullptr;
        const std::size_t scalar_count = weights.size() + lambda.size() + force.size();
        const std::size_t overlap_count = b.size() + db.size();
        ASSERT_EQ(cudaSuccess, cudaMallocManaged(&metadata, 6 * sizeof(int)));
        ASSERT_EQ(cudaSuccess, cudaMallocManaged(&scalars, scalar_count * sizeof(Real)));
        ASSERT_EQ(cudaSuccess, cudaMallocManaged(&overlaps, overlap_count * sizeof(std::complex<Real>)));
        std::copy(nh, nh + 2, metadata);
        std::copy(na, na + 2, metadata + 2);
        std::copy(isk, isk + 2, metadata + 4);
        Real* d_lambda = scalars + weights.size();
        Real* d_force = d_lambda + lambda.size();
        std::complex<Real>* d_db = overlaps + b.size();
        std::copy(weights.begin(), weights.end(), scalars);
        std::copy(lambda.begin(), lambda.end(), d_lambda);
        std::copy(force.begin(), force.end(), d_force);
        std::copy(b.begin(), b.end(), overlaps);
        std::copy(db.begin(), db.end(), d_db);
        const int* d_na = metadata + 2;
        // The launch wrapper reads isk on the host; use the host fixture as production does.
        hamilt::cal_force_nl_op<Real, base_device::DEVICE_GPU>()(
            nullptr, occupied, wg_nc, 2, 3, nbands, ik, nkb, npol, metadata, d_na, tpiba,
            scalars, d_lambda, isk, overlaps, d_db, d_force);
        ASSERT_EQ(cudaSuccess, cudaDeviceSynchronize());
        std::copy(d_force, d_force + force.size(), force.begin());
        EXPECT_EQ(cudaSuccess, cudaFree(overlaps));
        EXPECT_EQ(cudaSuccess, cudaFree(scalars));
        EXPECT_EQ(cudaSuccess, cudaFree(metadata));
    }
#endif
};

template <typename Real>
void check_cpu(const double tolerance)
{
    DeltaSpinFixture<Real> fixture;
    for (int axis = 0; axis < 3; ++axis)
    {
        fixture.prepare(2, axis);
        std::vector<Real> force(6, 0.125);
        fixture.cpu(2, force);
        fixture.check(force, 2, tolerance);
    }
    fixture.prepare(1, 2);
    for (int spin = 0; spin < 2; ++spin)
    {
        fixture.isk[1] = spin;
        std::vector<Real> force(6, 0.125);
        fixture.cpu(1, force);
        fixture.check(force, 1, tolerance);
    }
}

TEST(DeltaSpinForce, CpuDoubleEnergyDerivative)
{
    check_cpu<double>(1e-10);
}

TEST(DeltaSpinForce, CpuFloatEnergyDerivative)
{
    check_cpu<float>(2e-7);
}

#ifdef __CUDA
template <typename Real>
void check_gpu(const double tolerance)
{
    DeltaSpinFixture<Real> fixture;
    for (int npol : {1, 2})
    {
        for (int axis = 0; axis < 3; ++axis)
        {
            fixture.prepare(npol, axis);
            for (int spin = 0; spin < 2; ++spin)
            {
                fixture.isk[1] = spin;
                std::vector<Real> cpu(6, 0.125);
                std::vector<Real> gpu(6, 0.125);
                fixture.cpu(npol, cpu);
                fixture.gpu(npol, gpu);
                fixture.check(gpu, npol, tolerance);
                for (std::size_t i = 0; i < cpu.size(); ++i)
                {
                    EXPECT_NEAR(cpu[i], gpu[i], tolerance);
                }
            }
        }
    }
}

TEST(DeltaSpinForce, GpuDoubleEnergyDerivativeAndCpuParity)
{
    check_gpu<double>(1e-10);
}

TEST(DeltaSpinForce, GpuFloatEnergyDerivativeAndCpuParity)
{
    check_gpu<float>(2e-7);
}
#endif
} // namespace
