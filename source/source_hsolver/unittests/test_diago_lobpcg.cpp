#include "source_hsolver/diago_lobpcg.h"
#include "source_hsolver/rr_solver.h"
#include "source_base/module_container/base/third_party/lapack.h"
#include "source_base/module_device/device.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/parallel_cell.h"
#include "source_base/parallel_comm.h"
#include "source_base/parallel_device.h"
#include "source_basis/module_pw/test/test_tool.h"

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace
{
using Z = std::complex<double>;
using Locking = hsolver::LobpcgLocking;
#if defined(__CUDA) || defined(__ROCM)
using Device = base_device::DEVICE_GPU;
#else
using Device = base_device::DEVICE_CPU;
#endif
const ct::DeviceType device = ct::DeviceTypeToEnum<Device>::value;
class Lobpcg : public ::testing::TestWithParam<Locking> {};

hsolver::diag_comm_info comm_info()
{
    const auto world = ModuleBase::world_comm_domain();
#ifdef __MPI
    return {POOL_WORLD, world.rank(), world.size()};
#else
    return {0, 1};
#endif
}

template <typename T>
void sum(T* values, const int count)
{
#ifdef __MPI
    const auto comm = comm_info();
    Parallel_Common::reduce_data(values, count, comm.comm);
#endif
}

template <typename T>
void upload(T* destination, const std::vector<T>& source)
{
    base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>()(
        destination, source.data(), source.size());
}

template <typename T>
std::vector<T> download(const T* source, const int count)
{
    std::vector<T> result(count);
    base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, Device>()(
        result.data(), source, count);
    return result;
}

struct Options
{
    int n = 80;
    int bands = 6;
    int extra = 0;
    int npol = 2;
    int backend = 0;
    int exact = 0;
    int budget = 250;
    bool overlap = true;
    bool cluster = false;
    bool empty_rank = false;
};

// Independent dense H=L^H D L, S=L^H L; the same fixture runs on CPU/GPU.
// Only this test operator stages vectors on host. Production receives device pointers.
template <typename T>
class Problem : public hsolver::HSOperator<T, Device>
{
  public:
    using Real = typename GetTypeReal<T>::type;
    using Solver = hsolver::DiagoLOBPCG<T, Device>;
    Options options;
    hsolver::diag_comm_info comm;
    int first;
    int npw;
    int stride;
    int dim;
    std::vector<Z> h;
    std::vector<Z> s;
    std::vector<Z> reference_vectors;
    std::vector<double> reference_values;
    std::vector<T> input;
    std::vector<Real> prec;
    std::vector<Real> eigenvalues;
    ct::Tensor psi;
    hsolver::rr::ParallelRREigensolver parallel;
    hsolver::rr::GpuRREigensolver gpu;
    std::unique_ptr<Solver> solver;
    mutable int h_calls = 0;
    mutable int s_calls = 0;
    mutable int h_columns = 0;

    Problem(const Options& configuration, const Locking locking, const hsolver::diag_comm_info& domain)
        : options(configuration), comm(domain)
    {
        const int n = options.n;
        int local = n / comm.nproc;
        first = comm.rank * local;
        if (options.empty_rank && comm.nproc > 1)
        {
            local = comm.rank == 0 ? 0 : n / (comm.nproc - 1);
            first = std::max(0, comm.rank - 1) * local;
        }
        npw = local / options.npol;
        stride = npw + 3;
        dim = options.npol * stride;
        h.assign(n * n, Z(0));
        s.assign(n * n, Z(0));
        std::vector<Z> l(n * n, Z(0));
        for (int j = 0; j < n; ++j)
        {
            l[j + j * n] = 0.8 + 0.015 * j;
            for (int i = 0; i < j; ++i)
            {
                l[i + j * n] = Z(0.08 * std::sin(i + 2 * j), 0.06 * std::cos(2 * i + j));
            }
        }
        for (int j = 0; j < n; ++j)
        {
            for (int i = 0; i < n; ++i)
            {
                for (int k = 0; k < n; ++k)
                {
                    const double energy = options.cluster && k >= 3 && k < 8 ? 4.0 : k + 1.0;
                    const Z product = std::conj(l[k + i * n]) * l[k + j * n];
                    h[i + j * n] += energy * product;
                    s[i + j * n] += product;
                }
                if (!options.overlap) { s[i + j * n] = i == j ? 1 : 0; }
            }
        }
        reference_vectors = h;
        std::vector<Z> metric(s);
        reference_values.resize(n);
        std::vector<Z> work(2 * n);
        std::vector<double> real_work(3 * n);
        int info = 0;
        ct::lapackConnector::hegv(1, 'V', 'U', n, reference_vectors.data(), n, metric.data(), n,
                                 reference_values.data(), work.data(), 2 * n, real_work.data(), info);
        EXPECT_EQ(info, 0);
        input.assign(dim * options.bands, T(7, 9)); // Deliberately dirty input padding.
        for (int j = 0; j < options.bands; ++j)
        {
            for (int row = 0; row < dim; ++row)
            {
                if (!physical(row)) { continue; }
                const int global = global_row(row);
                Z value(std::sin(0.17 * (global + 1) * (j + 2)),
                        0.3 * std::cos(0.11 * (global + 3) * (j + 1)));
                if (j < options.exact) { value = reference_vectors[global + j * n]; }
                input[row + j * dim] = static_cast<T>(value);
            }
        }
        const auto type = ct::DataTypeToEnum<T>::value;
        psi = ct::Tensor(type, device, {dim, options.bands});
        prec.assign(dim, 1);
        eigenvalues.resize(options.bands);
        const int block = options.bands + options.extra;
        const hsolver::LobpcgBasisLayout layout = {npw, options.npol};
        solver.reset(new Solver(prec.data(), options.bands, dim, block, layout, comm,
                                options.backend, 2, options.overlap, locking, parallel, gpu));
        upload(psi.data<T>(), input);
    }

    bool physical(const int row) const { return row % stride < npw; }
    int global_row(const int row) const { return first + row / stride * npw + row % stride; }
    void update_k(const int) override {}
    void hpsi(const T* x, T* y, const int ld, const int count) const override
    {
        ++h_calls;
        h_columns += count;
        apply(h, x, y, ld, count);
    }
    void spsi(const T* x, T* y, const int ld, const int count) const override
    {
        ++s_calls;
        apply(s, x, y, ld, count);
    }
    void apply(const std::vector<Z>& matrix, const T* x, T* y, const int ld, const int count) const
    {
        const int elements = ld * count;
        const auto local = download(x, elements);
        std::vector<T> output(elements, T(91, -13)); // Operators need not clear output padding.
        for (int band = 0; band < count; ++band)
        {
            std::vector<Z> full(options.n, Z(0));
            for (int row = 0; row < dim; ++row)
            {
                if (physical(row)) { full[global_row(row)] = local[row + band * ld]; }
                else { EXPECT_EQ(local[row + band * ld], T(0)); }
            }
            sum(full.data(), options.n);
            for (int row = 0; row < dim; ++row)
            {
                if (!physical(row)) { continue; }
                Z value(0);
                for (int k = 0; k < options.n; ++k)
                {
                    value += matrix[global_row(row) + k * options.n] * full[k];
                }
                output[row + band * ld] = static_cast<T>(value);
            }
        }
        upload(y, output);
    }
    bool solve(const double tolerance, const std::vector<double>& thresholds)
    {
        return solver->diag(*this, psi.data<T>(), dim, eigenvalues.data(), tolerance, options.budget, thresholds);
    }
    // Check reference eigenvalues and the physical S Gram independently of solver projections.
    void check(const bool converged, const double eig_tol, const double gram_tol)
    {
        const auto x = download(psi.data<T>(), dim * options.bands);
        std::vector<Z> full(options.n * options.bands, Z(0));
        for (int j = 0; j < options.bands; ++j)
        {
            for (int row = 0; row < dim; ++row)
            {
                if (physical(row)) { full[global_row(row) + j * options.n] = x[row + j * dim]; }
                else { EXPECT_EQ(x[row + j * dim], T(0)); }
            }
        }
        const int count = full.size();
        sum(full.data(), count);
        std::vector<Z> sx(count, Z(0));
        for (int j = 0; j < options.bands; ++j)
        {
            for (int i = 0; i < options.n; ++i)
            {
                for (int k = 0; k < options.n; ++k)
                {
                    sx[i + j * options.n] += s[i + k * options.n] * full[k + j * options.n];
                }
            }
            EXPECT_TRUE(std::isfinite(eigenvalues[j]));
            if (converged)
            {
                EXPECT_NEAR(eigenvalues[j], reference_values[j], eig_tol);
            }
            for (int k = 0; k < options.bands; ++k)
            {
                Z gram(0);
                for (int i = 0; i < options.n; ++i)
                {
                    gram += std::conj(full[i + k * options.n]) * sx[i + j * options.n];
                }
                const Z expected = k == j ? Z(1) : Z(0);
                EXPECT_LT(std::abs(gram - expected), gram_tol);
            }
        }
    }
};

template <typename T>
void check_case(const Options& options, const Locking locking)
{
    Problem<T> p(options, locking, comm_info());
    const bool single = sizeof(T) == sizeof(std::complex<float>);
    const double tolerance = single ? 1e-7 : 1e-13;
    const bool converged = p.solve(tolerance, {});
    EXPECT_TRUE(converged);
    const double eig_tol = single ? 2e-4 : 1e-8;
    const double gram_tol = single ? 2e-5 : 1e-8;
    p.check(converged, eig_tol, gram_tol);
    if (options.overlap) { EXPECT_GT(p.s_calls, 0); }
    else { EXPECT_EQ(p.s_calls, 0); }
    const bool gpu_rr = device == ct::DeviceType::GpuDevice && hsolver::rr::GpuRREigensolver::available()
                        && p.comm.rank == 0 && options.backend == 0;
    if (gpu_rr) { EXPECT_GT(p.gpu.storage_bytes(), 0U); }
    else { EXPECT_EQ(p.gpu.storage_bytes(), 0U); }
}

template <typename T>
void check_eigenpair_identity(const Locking locking)
{
    Options o;
    o.n = 40;
    o.bands = 4;
    o.npol = 1;
    Problem<T> p(o, locking, comm_info());
    std::fill(p.h.begin(), p.h.end(), Z(0));
    std::fill(p.s.begin(), p.s.end(), Z(0));
    for (int i = 0; i < o.n; ++i)
    {
        const double metric = 1 + .05 * i;
        p.s[i + i * o.n] = metric;
        p.h[i + i * o.n] = (i + 1) * metric;
    }
    const bool single = sizeof(T) == sizeof(std::complex<float>);
    const double tolerance = single ? 1e-5 : 1e-11;
    ASSERT_TRUE(p.solve(tolerance, {}));
    EXPECT_GT(p.solver->iterations(), 0); // Exercise iterative output, not the exact-initial-state exit.
    const double eig_tol = single ? 2e-4 : 1e-8;
    const int elements = p.dim * o.bands;
    const auto x = download(p.psi.template data<T>(), elements);
    // Exact pairs are (j+1, e_j/sqrt(S_jj)), with unit gaps. Compare columns
    // up to phase to detect permutations and mixing of distinct eigenvectors.
    // Bounds are O(sqrt(eig_tol/gap)), allowing ordinary float/double roundoff.
    const double vector_tol = single ? 2e-2 : 2e-4;
    for (int j = 0; j < o.bands; ++j)
    {
        const double expected_eigenvalue = j + 1;
        EXPECT_NEAR(p.eigenvalues[j], expected_eigenvalue, eig_tol);
        double error_squared = 0;
        for (int row = 0; row < p.dim; ++row)
        {
            if (!p.physical(row)) { continue; }
            const int i = p.global_row(row);
            const double expected = i == j ? 1 : 0;
            const double amplitude = std::sqrt(1 + .05 * i) * std::abs(x[row + j * p.dim]);
            const double error = amplitude - expected;
            error_squared += error * error;
        }
        sum(&error_squared, 1);
        EXPECT_LT(std::sqrt(error_squared), vector_tol) << "band=" << j << ", single=" << single;
    }
}
} // namespace

namespace hsolver
{
// Only locking, projection-cache and recovery contracts require private access.
class LobpcgProjectionTest
{
  public:
    static void prefix(const Locking locking)
    {
        Problem<Z> p(Options{}, locking, comm_info());
        auto& solver = *p.solver;
        ct::Tensor previous(ct::DataType::DT_DOUBLE, ct::DeviceType::CpuDevice, {6});
        ct::Tensor threshold(ct::DataType::DT_DOUBLE, ct::DeviceType::CpuDevice, {6});
        for (int i = 0; i < 6; ++i)
        {
            solver.eig_.data<double>()[i] = previous.data<double>()[i] = i + 1;
            threshold.data<double>()[i] = 1e-9;
        }
        previous.data<double>()[5] += 1;
        int locked = 0;
        EXPECT_FALSE(solver.lock_converged(previous, threshold, locked));
        EXPECT_EQ(locked, 5);
        const int changed = locking == Locking::Soft ? 1 : 5;
        solver.eig_.data<double>()[changed] += 0.1;
        EXPECT_FALSE(solver.lock_converged(previous, threshold, locked));
        EXPECT_EQ(locked, locking == Locking::Soft ? 1 : 5);
        EXPECT_TRUE(solver.lock_converged(previous, threshold, locked));
        EXPECT_EQ(locked, 6);
    }

    static void cache(const int backend)
    {
        Options o;
        o.extra = 2;
        o.backend = backend;
        Problem<Z> p(o, Locking::Hard, comm_info());
        auto& solver = *p.solver;
        const int m = 9;
        const int q = 6;
        const auto basis = solver.basis_images(0);
        std::vector<Z> x(p.dim * m, Z(0));
        for (int j = 0; j < m; ++j)
        {
            for (int i = 0; i < p.dim; ++i)
            {
                if (p.physical(i))
                {
                    const int row = p.global_row(i);
                    x[i + j * p.dim] = Z(std::sin(.13 * (row + 1) * (j + 2)), std::cos(.17 * (row + 2) * (j + 1))) / 8.0;
                }
            }
        }
        upload(basis.x, x);
        p.hpsi(basis.x, basis.hx, p.dim, m);
        p.spsi(basis.x, basis.sx, p.dim, m);
        solver.project_padding(basis.hx, m);
        solver.project_padding(basis.sx, m);
        ct::Tensor k;
        ct::Tensor g;
        solver.form_projected_matrices(basis, m, nullptr, k, g);
        const double nan = std::numeric_limits<double>::quiet_NaN();
        std::fill_n(solver.ritz_coefficients_.data<Z>(), solver.ritz_coefficients_.NumElements(), Z(nan));
        std::fill_n(solver.e_red_.data<double>(), m, nan);
        rr::SerialGvWorkspace work;
        ASSERT_TRUE(solver.solve_projected(k, g, m, q, work));
        for (int i = q; i < m; ++i) { EXPECT_TRUE(std::isnan(solver.e_red_.data<double>()[i])); }
        for (int i = m * q; i < solver.ritz_coefficients_.NumElements(); ++i)
        {
            EXPECT_TRUE(std::isnan(solver.ritz_coefficients_.data<Z>()[i].real()));
        }
        ct::Tensor cx(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, {m, q});
        std::copy_n(solver.ritz_coefficients_.data<Z>(), m * q, cx.data<Z>());
        const auto hx = download(basis.hx, p.dim * m);
        const auto sx = download(basis.sx, p.dim * m);
        for (const int first : {0, 2})
        {
            SCOPED_TRACE(first);
            upload(basis.x, x);
            upload(basis.hx, hx);
            upload(basis.sx, sx);
            const int active = q - first;
            const int np = solver.build_direction(basis, m, q, active, cx, g);
            ASSERT_GT(np, 0);
            const auto history = solver.basis_images(solver.n_max_);
            const auto actual_p = download(history.x, p.dim * np);
            const auto actual_hp = download(history.hx, p.dim * np);
            const auto actual_sp = download(history.sx, p.dim * np);
            const Z* cp = solver.u_p_.data<Z>();
            for (int j = 0; j < np; ++j)
            {
                double norm = 0;
                for (int row = 0; row < p.dim; ++row)
                {
                    Z expected(0);
                    Z expected_h(0);
                    Z expected_s(0);
                    for (int l = 0; l < m; ++l)
                    {
                        expected += x[row + l * p.dim] * cp[l + j * m];
                        expected_h += hx[row + l * p.dim] * cp[l + j * m];
                        expected_s += sx[row + l * p.dim] * cp[l + j * m];
                    }
                    EXPECT_LT(std::abs(actual_p[row + j * p.dim] - expected), 1e-11);
                    EXPECT_LT(std::abs(actual_hp[row + j * p.dim] - expected_h), 1e-11);
                    EXPECT_LT(std::abs(actual_sp[row + j * p.dim] - expected_s), 1e-11);
                    norm += std::norm(actual_p[row + j * p.dim]);
                }
                sum(&norm, 1);
                EXPECT_NEAR(norm, 1.0, 1e-12);
                if (first == 0)
                {
                    for (int i = 0; i < q; ++i) { EXPECT_EQ(cp[i + j * m], Z(0)); }
                }
                for (int i = 0; i < first && p.comm.rank == 0; ++i)
                {
                    Z overlap(0);
                    for (int a = 0; a < m; ++a)
                    {
                        for (int b = 0; b < m; ++b)
                        {
                            overlap += std::conj(cx.data<Z>()[a + i * m]) * g.data<Z>()[a + b * m] * cp[b + j * m];
                        }
                    }
                    EXPECT_LT(std::abs(overlap), 1e-10);
                }
            }
            typename Problem<Z>::Solver::ProjectionCache retained;
            const double* values = solver.e_red_.data<double>() + first;
            const int zero_prefix = first == 0 ? q : 0;
            solver.retain_pair(k, g, m, cx, first, active, np, values, zero_prefix, retained);
            const int keep = active + np;
            std::vector<Z> next(p.dim * (keep + 1), Z(0));
            for (int j = 0; j < active; ++j)
            {
                for (int row = 0; row < p.dim; ++row)
                {
                    for (int l = 0; l < m; ++l)
                    {
                        next[row + j * p.dim] += x[row + l * p.dim] * cx.data<Z>()[l + (first + j) * m];
                    }
                }
            }
            std::copy(actual_p.begin(), actual_p.end(), next.begin() + active * p.dim);
            std::copy_n(x.data(), p.dim, next.data() + keep * p.dim);
            upload(basis.x, next);
            const int count = keep + 1;
            p.hpsi(basis.x, basis.hx, p.dim, count);
            p.spsi(basis.x, basis.sx, p.dim, count);
            solver.project_padding(basis.hx, count);
            solver.project_padding(basis.sx, count);
            for (const int columns : {keep, count})
            {
                ct::Tensor cached_k;
                ct::Tensor cached_g;
                ct::Tensor full_k;
                ct::Tensor full_g;
                solver.form_projected_matrices(basis, columns, &retained, cached_k, cached_g);
                solver.form_projected_matrices(basis, columns, nullptr, full_k, full_g);
                for (int i = 0; i < full_k.NumElements(); ++i)
                {
                    EXPECT_LT(std::abs(cached_k.data<Z>()[i] - full_k.data<Z>()[i]), 2e-12);
                    EXPECT_LT(std::abs(cached_g.data<Z>()[i] - full_g.data<Z>()[i]), 2e-12);
                }
            }
            solver.retain_pair(k, g, m, cx, first, active, 0, values, 0, retained);
            EXPECT_EQ(retained.columns, active);
            for (int j = 0; j < active && p.comm.rank == 0; ++j)
            {
                for (int i = 0; i < active; ++i)
                {
                    const Z metric = i == j ? 1 : 0;
                    const Z energy = i == j ? values[i] : 0;
                    EXPECT_LT(std::abs(retained.g.data<Z>()[i + j * active] - metric), 2e-12);
                    EXPECT_LT(std::abs(retained.k.data<Z>()[i + j * active] - energy), 2e-12);
                }
            }
        }
    }

    static void recovery(const bool dependent)
    {
        Options o;
        o.bands = 4;
        Problem<Z> p(o, Locking::Hard, comm_info());
        auto& solver = *p.solver;
        std::fill(p.h.begin(), p.h.end(), Z(0));
        std::fill(p.s.begin(), p.s.end(), Z(0));
        for (int i = 0; i < o.n; ++i)
        {
            p.s[i + i * o.n] = 1 + .01 * i;
            p.h[i + i * o.n] = (i + 1.0) * p.s[i + i * o.n];
        }
        const auto basis = solver.basis_images(0);
        std::vector<Z> x(p.dim * 7, Z(0));
        for (int row = 0; row < p.dim; ++row)
        {
            if (!p.physical(row)) { continue; }
            const int i = p.global_row(row);
            for (int j = i; j < 4; ++j)
            {
                x[row + j * p.dim] = i == j ? Z(.5 + .1 * i) : Z(.08, .03);
            } // Full-rank reference with nonunit diagonal and nonzero off-diagonal S overlaps.
            for (int j = 4; j < 7; ++j)
            {
                const int frequency = dependent && j == 6 ? 4 : j;
                x[row + j * p.dim] = Z(std::sin(.13 * (i + 1) * frequency), std::cos(.19 * i * frequency));
            }
        }
        upload(basis.x, x);
        p.hpsi(basis.x, basis.hx, p.dim, 7);
        p.spsi(basis.x, basis.sx, p.dim, 7);
        solver.project_padding(basis.hx, 7);
        solver.project_padding(basis.sx, 7);
        ct::Tensor k;
        ct::Tensor g;
        solver.form_projected_matrices(basis, 7, nullptr, k, g);
        if (p.comm.rank == 0) { g.data<Z>()[0] = -1; }
        rr::SerialGvWorkspace work;
        EXPECT_FALSE(solver.solve_projected(k, g, 7, 4, work));
        EXPECT_NE(solver.error_message().find("dimension=7"), std::string::npos);
        EXPECT_NE(solver.error_message().find("wanted=4"), std::string::npos);
        int columns = 7;
        ASSERT_TRUE(solver.recover_projected(p, basis, columns, 4, k, g, work));
        EXPECT_TRUE(solver.error_message().empty());
        // Exercise rank repair independently; the preceding recovery may have compacted the basis.
        upload(basis.x, x);
        p.hpsi(basis.x, basis.hx, p.dim, 7);
        solver.project_padding(basis.hx, 7);
        const int kept = solver.repair_basis(p, 3);
        const int expected_rank = dependent ? 2 : 3;
        ASSERT_EQ(kept, expected_rank);
        const int count = 4 + kept;
        const auto repaired = download(basis.x, p.dim * count);
        const auto repaired_s = download(basis.sx, p.dim * count);
        const auto repaired_h = download(basis.hx, p.dim * count);
        for (int i = 0; i < p.dim * 4; ++i) { EXPECT_EQ(repaired[i], x[i]); }
        for (int j = 0; j < count; ++j)
        {
            for (int row = 0; row < p.dim; ++row)
            {
                if (!p.physical(row)) { continue; }
                const int i = p.global_row(row);
                EXPECT_LT(std::abs(repaired_s[row + j * p.dim] - (1 + .01 * i) * repaired[row + j * p.dim]), 1e-10);
                EXPECT_LT(std::abs(repaired_h[row + j * p.dim] - (i + 1.0) * repaired_s[row + j * p.dim]), 1e-10);
            }
            for (int i = 0; i < count; ++i)
            {
                Z gram(0);
                for (int row = 0; row < p.dim; ++row)
                {
                    if (p.physical(row)) { gram += std::conj(repaired[row + i * p.dim]) * repaired_s[row + j * p.dim]; }
                }
                sum(&gram, 1);
                Z expected = i == j && i >= 4 ? 1 : 0;
                if (i < 4 && j < 4)
                {
                    for (int row = 0; row < p.dim; ++row)
                    {
                        const double metric = 1 + .01 * p.global_row(row);
                        if (p.physical(row)) { expected += std::conj(x[row + i * p.dim]) * metric * x[row + j * p.dim]; }
                    }
                    sum(&expected, 1);
                }
                EXPECT_LT(std::abs(gram - expected), 1e-10);
            }
        }
        // Invalid reference X must fail without modifying its columns.
        std::fill_n(x.data(), p.dim, Z(0));
        upload(basis.x, x);
        EXPECT_LT(solver.repair_basis(p, 3), 0);
        EXPECT_EQ(download(basis.x, p.dim * 4), std::vector<Z>(x.begin(), x.begin() + p.dim * 4));
    }

    static void fatal(const std::string& scenario, const Locking locking)
    {
        Options o;
        o.n = 120;
        o.bands = 4;
        o.npol = 1;
        o.extra = scenario == "block" ? -1 : 0;
        o.backend = scenario == "unsupported_backend" ? 1 : 0;
        o.overlap = scenario == "overlap_negative" || scenario == "overlap_zero";
        if (scenario == "reused_invalid")
        {
            o.n = 80;
            o.bands = 6;
            o.extra = 5;
            o.npol = 2;
        }
        Problem<Z> p(o, locking, comm_info());
        if (scenario == "dimension")
        {
            const LobpcgBasisLayout tiny = {1, 1};
            DiagoLOBPCG<Z, Device> solver(p.prec.data(), 4, p.dim, 4, tiny, p.comm, 0, 2, false, locking, p.parallel, p.gpu);
            solver.diag(p, p.psi.data<Z>(), p.dim, p.eigenvalues.data(), 1e-8, 50, {});
        }
        else if (scenario == "search")
        {
            const double nan = std::numeric_limits<double>::quiet_NaN();
            std::fill(p.input.begin(), p.input.end(), Z(nan));
            upload(p.psi.data<Z>(), p.input);
            const DiagoLOBPCG<Z, Device>::VectorImages images = {p.psi.data<Z>(), nullptr, nullptr};
            p.solver->normalize_columns(images, o.bands, nullptr, 0);
        }
        else
        {
            if (scenario == "reused_invalid") { ASSERT_TRUE(p.solve(1e-8, {})); }
            if (o.overlap)
            {
                std::fill(p.s.begin(), p.s.end(), Z(0));
                if (scenario == "overlap_negative")
                {
                    for (int i = 0; i < o.n; ++i) { p.s[i + i * o.n] = -1; }
                }
            }
            else if (scenario == "dependent" || scenario == "reused_invalid")
            {
                std::fill(p.input.begin(), p.input.end(), Z(.3, .1));
                upload(p.psi.data<Z>(), p.input);
            }
            p.solve(1e-8, {});
        }
    }
};
} // namespace hsolver

namespace
{
TEST_P(Lobpcg, StandardGeneralizedSpinorAndBuffer)
{
    for (const bool overlap : {false, true})
    {
        for (const int extra : {0, 5})
        {
            SCOPED_TRACE(::testing::Message() << "S=" << overlap << " extra=" << extra);
            Options o;
            o.overlap = overlap;
            o.npol = extra == 0 ? 1 : 2;
            o.extra = extra;
            check_case<Z>(o, GetParam());
        }
    }
}
TEST_P(Lobpcg, ClusterAndPartiallyExact)
{
    for (const int exact : {0, 2})
    {
        Options o;
        o.cluster = true;
        o.extra = exact == 0 ? 4 : 0;
        o.exact = exact;
        check_case<Z>(o, GetParam());
    }
}
TEST_P(Lobpcg, SinglePrecision)
{
    for (const bool cluster : {false, true})
    {
        Options o;
        o.cluster = cluster;
        o.extra = cluster ? 2 : 0;
        check_case<std::complex<float>>(o, GetParam());
    }
}
TEST_P(Lobpcg, EigenvectorsMatchEigenvalues)
{
    check_eigenpair_identity<Z>(GetParam());
    check_eigenpair_identity<std::complex<float>>(GetParam());
}
TEST_P(Lobpcg, DistributedRitzAndEmptyRank)
{
    if (comm_info().nproc == 1) { GTEST_SKIP() << "Requires multiple ranks"; }
    for (const bool empty : {false, true})
    {
        Options o;
        o.backend = 2;
        o.empty_rank = empty;
        if (empty) { o.n = 24 * (comm_info().nproc - 1); }
        check_case<Z>(o, GetParam());
        if (!empty) { check_case<std::complex<float>>(o, GetParam()); }
    }
}
TEST_P(Lobpcg, BudgetAndReuse)
{
    Options o;
    o.budget = 1;
    Problem<Z> p(o, GetParam(), comm_info());
    EXPECT_FALSE(p.solve(1e-13, {}));
    EXPECT_EQ(p.solver->iterations(), 1);
    p.check(false, 1e-8, 1e-8);
    p.options.budget = 250;
    EXPECT_TRUE(p.solve(1e-13, {}));
    p.check(true, 1e-8, 1e-8);
}
TEST(LobpcgInitialization, ExactNonorthogonalInitialization)
{
    Options o;
    o.bands = 4;
    o.exact = 4;
    Problem<Z> p(o, Locking::Hard, comm_info());
    const auto original = p.input;
    for (int j = 0; j < o.bands; ++j)
    {
        for (int row = 0; row < p.dim; ++row)
        {
            if (!p.physical(row)) { continue; }
            Z value(0);
            for (int k = 0; k <= j; ++k)
            {
                const Z mix = k == j ? Z(1 + .3 * j, .1 * j) : Z(.2, -.1);
                value += original[row + k * p.dim] * mix;
            }
            p.input[row + j * p.dim] = value;
        }
    }
    upload(p.psi.data<Z>(), p.input);
    EXPECT_TRUE(p.solve(1e-10, {}));
    EXPECT_LE(p.solver->iterations(), 1);
    EXPECT_EQ(p.h_calls, 1 + p.solver->iterations());
    EXPECT_EQ(p.s_calls, 1 + p.solver->iterations());
    p.check(true, 1e-10, 1e-10);
}
TEST_P(Lobpcg, IdentityMetricAgreesWithImplicitOverlap)
{
    Options o;
    o.overlap = false;
    Problem<Z> implicit(o, GetParam(), comm_info());
    ASSERT_TRUE(implicit.solve(1e-12, {}));
    EXPECT_EQ(implicit.s_calls, 0);
    o.overlap = true;
    Problem<Z> explicit_s(o, GetParam(), comm_info());
    explicit_s.s = implicit.s;
    explicit_s.reference_values = implicit.reference_values;
    ASSERT_TRUE(explicit_s.solve(1e-12, {}));
    EXPECT_GT(explicit_s.s_calls, 0);
    explicit_s.check(true, 1e-8, 1e-8);
    for (int i = 0; i < o.bands; ++i)
    {
        EXPECT_NEAR(explicit_s.eigenvalues[i], implicit.eigenvalues[i], 1e-9);
    }
}
TEST(LobpcgInitialization, ExactlyRepresentableState)
{
    Options o;
    o.bands = 4;
    Problem<Z> p(o, Locking::Hard, comm_info());
    std::fill(p.h.begin(), p.h.end(), Z(0));
    std::fill(p.s.begin(), p.s.end(), Z(0));
    std::fill(p.input.begin(), p.input.end(), Z(0));
    for (int i = 0; i < o.n; ++i)
    {
        const double metric = std::pow(4.0, i % 5 - 2);
        p.s[i + i * o.n] = metric;
        p.h[i + i * o.n] = (i + 1) * metric;
        p.reference_values[i] = i + 1;
    }
    for (int row = 0; row < p.dim; ++row)
    {
        if (!p.physical(row)) { continue; }
        const int i = p.global_row(row);
        if (i < o.bands) { p.input[row + i * p.dim] = 1 / std::sqrt(p.s[i + i * o.n].real()); }
    }
    upload(p.psi.data<Z>(), p.input);
    EXPECT_TRUE(p.solve(1e-12, {}));
    EXPECT_EQ(p.solver->iterations(), 0);
    EXPECT_EQ(p.h_calls, 1);
    p.check(true, 1e-12, 1e-12);
}
TEST_P(Lobpcg, BandThresholdOverrideAndRepeatability)
{
    Options o;
    o.extra = 4;
    Problem<Z> p(o, GetParam(), comm_info());
    std::vector<Z> vectors;
    std::vector<double> eigenvalues;
    int iterations = 0;
    int h_columns = 0;
    for (int variant = 0; variant < 3; ++variant)
    {
        double scalar = 1e-13;
        std::vector<double> thresholds;
        if (variant == 1) { scalar = 1e-3; thresholds.assign(o.bands, 1e-13); }
        if (variant == 2) { thresholds.assign(2, 1e-13); }
        upload(p.psi.data<Z>(), p.input);
        p.h_columns = 0;
        ASSERT_TRUE(p.solve(scalar, thresholds));
        const auto actual = download(p.psi.data<Z>(), p.dim * o.bands);
        if (variant == 0)
        {
            vectors = actual;
            eigenvalues = p.eigenvalues;
            iterations = p.solver->iterations();
            h_columns = p.h_columns;
        }
        else
        {
            EXPECT_EQ(p.solver->iterations(), iterations);
            EXPECT_EQ(p.h_columns, h_columns);
            for (int i = 0; i < o.bands; ++i) { EXPECT_NEAR(p.eigenvalues[i], eigenvalues[i], 1e-12); }
            for (int i = 0; i < actual.size(); ++i) { EXPECT_LT(std::abs(actual[i] - vectors[i]), 1e-10); }
        }
    }
}
TEST_P(Lobpcg, PrefixReactivationPolicy) { hsolver::LobpcgProjectionTest::prefix(GetParam()); }
TEST(LobpcgSubspace, RetainAndPartialRitz)
{
    hsolver::LobpcgProjectionTest::cache(0);
    if (comm_info().nproc > 1) { hsolver::LobpcgProjectionTest::cache(2); }
}
TEST(LobpcgSubspace, MetricFailureAndRankRepair)
{
    hsolver::LobpcgProjectionTest::recovery(false);
    hsolver::LobpcgProjectionTest::recovery(true);
}
} // namespace

int main(int argc, char** argv)
{
#ifdef __MPI
    int nproc = 1;
    int rank = 0;
    int in_pool = 0;
    int kpar = 1;
    int pool = 0;
    int pool_rank = 0;
    setupmpi(argc, argv, nproc, rank);
    divide_pools(nproc, rank, in_pool, kpar, pool, pool_rank);
#endif
#if defined(__CUDA) || defined(__ROCM)
    if (!base_device::information::probe_gpu_availability())
    {
#ifdef __MPI
        finishmpi();
#endif
        return 77;
    }
    ModuleBase::createGpuBlasHandle();
    ct::kernels::createGpuBlasHandle();
#endif
    if (argc == 4 && std::string(argv[1]) == "--lobpcg-fatal")
    {
        const Locking locking = std::string(argv[3]) == "soft" ? Locking::Soft : Locking::Hard;
        hsolver::LobpcgProjectionTest::fatal(argv[2], locking);
        std::cerr << "Fatal LOBPCG test unexpectedly returned" << std::endl;
        return 99;
    }
    testing::InitGoogleTest(&argc, argv);
    auto& listeners = testing::UnitTest::GetInstance()->listeners();
    if (comm_info().rank != 0) { delete listeners.Release(listeners.default_result_printer()); }
    const int result = RUN_ALL_TESTS();
#if defined(__CUDA) || defined(__ROCM)
    ct::kernels::destroyGpuBlasHandle();
    ModuleBase::destoryBLAShandle();
#endif
#ifdef __MPI
    finishmpi();
#endif
    return result;
}

INSTANTIATE_TEST_SUITE_P(LockingModes, Lobpcg, ::testing::Values(Locking::Hard, Locking::Soft));
