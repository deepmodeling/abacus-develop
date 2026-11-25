#include "source_hsolver/diago_lobpcg.h"
#include "source_hsolver/dense_eigensolver.h"
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
#include <memory>
#include <type_traits>
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
    const ModuleBase::CommunicationDomain world = ModuleBase::world_comm_domain();
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
    const hsolver::diag_comm_info comm = comm_info();
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
    int npol = 1;
    int backend = 0;
    bool overlap = true;
};

// Small H/S matrices applied on host, with distributed rows and CPU/GPU vector storage.
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
    hsolver::dense::ParallelDenseEigensolver parallel;
    hsolver::dense::GpuDenseEigensolver gpu;
    std::unique_ptr<Solver> solver;
    mutable int h_calls = 0;
    mutable int s_calls = 0;

    Problem(const Options& configuration, const Locking locking, const hsolver::diag_comm_info& domain)
        : options(configuration), comm(domain)
    {
        const int n = options.n;
        int local = n / comm.nproc;
        first = comm.rank * local;
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
                    const double energy = k + 1.0;
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
        const int lwork = 2 * n;
        std::vector<Z> work(lwork);
        std::vector<double> real_work(3 * n);
        int info = 0;
        ct::lapackConnector::hegv(1, 'V', 'U', n, reference_vectors.data(), n, metric.data(), n,
                                 reference_values.data(), work.data(), lwork, real_work.data(), info);
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
                input[row + j * dim] = static_cast<T>(value);
            }
        }
        const ct::DataType type = ct::DataTypeToEnum<T>::value;
        psi = ct::Tensor(type, device, {dim, options.bands});
        prec.assign(dim, 1);
        eigenvalues.resize(options.bands);
        const int block = options.bands + options.extra;
        const hsolver::LobpcgBasisLayout layout = {npw, options.npol};
        solver.reset(new Solver(prec.data(), options.bands, dim, block, layout, comm,
                                options.backend, 2, options.overlap, locking, parallel, gpu));
        upload(psi.data<T>(), input);
    }

    void diagonal_problem()
    {
        const int n = options.n;
        std::fill(h.begin(), h.end(), Z(0));
        std::fill(s.begin(), s.end(), Z(0));
        std::fill(reference_vectors.begin(), reference_vectors.end(), Z(0));
        for (int i = 0; i < n; ++i)
        {
            const double metric = options.overlap ? 1 + .05 * i : 1;
            s[i + i * n] = metric;
            h[i + i * n] = (i + 1) * metric;
            reference_values[i] = i + 1;
            reference_vectors[i + i * n] = 1 / std::sqrt(metric);
        }
    }
    void exact_initial(const int bands)
    {
        for (int j = 0; j < bands; ++j)
        {
            for (int row = 0; row < dim; ++row)
            {
                if (physical(row)) { input[row + j * dim] = T(reference_vectors[global_row(row) + j * options.n]); }
            }
        }
        upload(psi.data<T>(), input);
    }
    bool physical(const int row) const { return row % stride < npw; }
    int global_row(const int row) const { return first + row / stride * npw + row % stride; }
    void update_k(const int) override {}
    void hpsi(const T* x, T* y, const int ld, const int count) const override
    {
        ++h_calls;
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
        const std::vector<T> local = download(x, elements);
        std::vector<T> output(elements, T(91, -13)); // Sentinel padding outside physical rows.
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
    bool solve(const double tolerance, const int budget, const std::vector<double>& thresholds)
    {
        return solver->diag(*this, psi.data<T>(), dim, eigenvalues.data(), tolerance, budget, thresholds);
    }
    // Check reference eigenvalues and the physical S Gram independently of solver projections.
    void check(const bool converged, const double eig_tol, const double gram_tol)
    {
        const std::vector<T> x = download(psi.data<T>(), dim * options.bands);
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
    const bool single = std::is_same<T, std::complex<float>>::value;
    const double tolerance = single ? 1e-7 : 1e-13;
    const bool converged = p.solve(tolerance, 250, {});
    EXPECT_TRUE(converged);
    const double eig_tol = single ? 2e-4 : 1e-8;
    const double gram_tol = single ? 2e-5 : 1e-8;
    p.check(converged, eig_tol, gram_tol);
    if (options.overlap) { EXPECT_GT(p.s_calls, 0); }
    else { EXPECT_EQ(p.s_calls, 0); }
    const bool gpu_rr = device == ct::DeviceType::GpuDevice && hsolver::dense::GpuDenseEigensolver::available()
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
    p.diagonal_problem();
    const bool single = std::is_same<T, std::complex<float>>::value;
    const double tolerance = single ? 1e-5 : 1e-11;
    ASSERT_TRUE(p.solve(tolerance, 250, {}));
    EXPECT_GT(p.solver->iterations(), 0); // Exercise the iterative solve path.
    const double eig_tol = single ? 2e-4 : 1e-8;
    const int elements = p.dim * o.bands;
    const std::vector<T> x = download(p.psi.template data<T>(), elements);
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
template <typename T>
void check_exact_initial(const Locking locking)
{
    Options o;
    o.bands = 4;
    o.overlap = false;
    Problem<T> p(o, locking, comm_info());
    p.diagonal_problem();
    p.exact_initial(o.bands);
    ASSERT_TRUE(p.solve(1e-5, 50, {}));
    EXPECT_EQ(p.solver->iterations(), 0);
    p.check(true, 1e-6, 1e-6);
}
} // namespace

namespace hsolver
{
class LobpcgProjectionTest
{
  public:
    static void prefix(const Locking locking)
    {
        Problem<Z> p(Options{}, locking, comm_info());
        Problem<Z>::Solver& solver = *p.solver;
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

    // Columns e_i/sqrt(S_ii) give G=I and K=diag(1,...,columns).
    static std::vector<Z> coordinate_basis(Problem<Z>& p, const int columns)
    {
        p.diagonal_problem();
        std::vector<Z> x(p.dim * columns, Z(0));
        for (int row = 0; row < p.dim; ++row)
        {
            if (!p.physical(row)) { continue; }
            const int i = p.global_row(row);
            if (i < columns) { x[row + i * p.dim] = p.reference_vectors[i + i * p.options.n]; }
        }
        const Problem<Z>::Solver::VectorImages basis = p.solver->basis_images(0);
        upload(basis.x, x);
        p.hpsi(basis.x, basis.hx, p.dim, columns);
        p.spsi(basis.x, basis.sx, p.dim, columns);
        p.solver->project_padding(basis.hx, columns);
        p.solver->project_padding(basis.sx, columns);
        return x;
    }

    static void direction(const int locked)
    {
        Problem<Z> p(Options{}, Locking::Hard, comm_info());
        Problem<Z>::Solver& solver = *p.solver;
        const std::vector<Z> x = coordinate_basis(p, 4);
        const Problem<Z>::Solver::VectorImages basis = solver.basis_images(0);
        const std::vector<Z> hx = download(basis.hx, p.dim * 4);
        const std::vector<Z> sx = download(basis.sx, p.dim * 4);
        ct::Tensor k;
        ct::Tensor g;
        solver.form_projected_matrices(basis, 4, nullptr, k, g);
        ct::Tensor c(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, {4, 2});
        c.zero();
        c.data<Z>()[0] = c.data<Z>()[2] = std::sqrt(.5);
        c.data<Z>()[4] = .5;
        c.data<Z>()[6] = -.5;
        c.data<Z>()[7] = std::sqrt(.5);
        const int active = 2 - locked;
        ASSERT_EQ(solver.build_direction(basis, 4, 2, active, c, g), active);
        const Problem<Z>::Solver::VectorImages history = solver.basis_images(solver.n_max_);
        const Z* inputs[] = {x.data(), hx.data(), sx.data()};
        const Z* outputs[] = {history.x, history.hx, history.sx};
        const Z* cp = solver.u_p_.data<Z>();
        for (int image = 0; image < 3; ++image)
        {
            const std::vector<Z> actual = download(outputs[image], p.dim * active);
            for (int j = 0; j < active; ++j)
            {
                for (int row = 0; row < p.dim; ++row)
                {
                    Z expected(0);
                    for (int i = 0; i < 4; ++i) { expected += inputs[image][row + i * p.dim] * cp[i + j * 4]; }
                    EXPECT_LT(std::abs(actual[row + j * p.dim] - expected), 1e-12);
                }
                // G=I: the projected direction is orthogonal to the locked Ritz vector.
                if (locked) { EXPECT_LT(std::abs(cp[j * 4] + cp[2 + j * 4]), 1e-12); }
                else { EXPECT_EQ(cp[j * 4], Z(0)); EXPECT_EQ(cp[1 + j * 4], Z(0)); }
            }
        }
    }

    static void retained_projection(const int zero_prefix)
    {
        Problem<Z> p(Options{}, Locking::Hard, comm_info());
        Problem<Z>::Solver& solver = *p.solver;
        const std::vector<Z> original = coordinate_basis(p, 5);
        const Problem<Z>::Solver::VectorImages basis = solver.basis_images(0);
        ct::Tensor k;
        ct::Tensor g;
        solver.form_projected_matrices(basis, 5, nullptr, k, g);
        ct::Tensor cx(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::CpuDevice, {5, 2});
        cx.zero();
        cx.data<Z>()[0] = cx.data<Z>()[6] = 1;
        Z* cp = solver.u_p_.data<Z>();
        std::fill_n(cp, 5, Z(0));
        const int first = std::max(1, zero_prefix);
        cp[first] = .6;
        cp[first + 1] = Z(0, .8);
        const double values[] = {1, 2};
        Problem<Z>::Solver::ProjectionCache retained;
        solver.retain_pair(k, g, 5, cx, 0, 2, 1, values, zero_prefix, retained);
        EXPECT_EQ(retained.columns, 3);
        // Form the physical [X,P,W] independently of the cached small matrices.
        std::vector<Z> next(p.dim * 4, Z(0));
        std::copy_n(original.data(), p.dim * 2, next.data());
        for (int row = 0; row < p.dim; ++row)
        {
            for (int i = 0; i < 5; ++i) { next[row + 2 * p.dim] += original[row + i * p.dim] * cp[i]; }
            next[row + 3 * p.dim] = original[row + 4 * p.dim];
        }
        upload(basis.x, next);
        p.hpsi(basis.x, basis.hx, p.dim, 4);
        p.spsi(basis.x, basis.sx, p.dim, 4);
        solver.project_padding(basis.hx, 4);
        solver.project_padding(basis.sx, 4);
        ct::Tensor cached_k;
        ct::Tensor cached_g;
        ct::Tensor full_k;
        ct::Tensor full_g;
        solver.form_projected_matrices(basis, 4, &retained, cached_k, cached_g);
        solver.form_projected_matrices(basis, 4, nullptr, full_k, full_g);
        for (int i = 0; i < full_k.NumElements(); ++i)
        {
            EXPECT_LT(std::abs(cached_k.data<Z>()[i] - full_k.data<Z>()[i]), 1e-12);
            EXPECT_LT(std::abs(cached_g.data<Z>()[i] - full_g.data<Z>()[i]), 1e-12);
        }
    }

    static void recovery()
    {
        Options o;
        o.bands = 4;
        Problem<Z> p(o, Locking::Hard, comm_info());
        coordinate_basis(p, 7);
        Problem<Z>::Solver& solver = *p.solver;
        const Problem<Z>::Solver::VectorImages basis = solver.basis_images(0);
        ct::Tensor k;
        ct::Tensor g;
        solver.form_projected_matrices(basis, 7, nullptr, k, g);
        if (p.comm.rank == 0) { g.data<Z>()[0] = -1; }
        dense::SerialGvWorkspace work;
        EXPECT_FALSE(solver.solve_projected(k, g, 7, 4, work));
        int columns = 7;
        ASSERT_TRUE(solver.recover_projected(p, basis, columns, 4, k, g, work));
        EXPECT_EQ(columns, 7);
        EXPECT_TRUE(solver.error_message_.empty());
        for (int i = 0; i < 4; ++i) { EXPECT_NEAR(solver.e_red_.data<double>()[i], i + 1, 1e-12); }
    }

    static void rank_repair()
    {
        Options o;
        o.bands = 4;
        Problem<Z> p(o, Locking::Hard, comm_info());
        std::vector<Z> x = coordinate_basis(p, 7);
        Problem<Z>::Solver& solver = *p.solver;
        const Problem<Z>::Solver::VectorImages basis = solver.basis_images(0);
        for (int row = 0; row < p.dim; ++row)
        {
            // Nonorthogonal, full-rank X and two independent increments plus a duplicate.
            x[row + p.dim] += .3 * x[row];
            x[row + 4 * p.dim] += .2 * x[row];
            x[row + 6 * p.dim] = x[row + 4 * p.dim];
        }
        upload(basis.x, x);
        p.hpsi(basis.x, basis.hx, p.dim, 7);
        solver.project_padding(basis.hx, 7);
        ASSERT_EQ(solver.repair_basis(p, 3), 2);
        const std::vector<Z> repaired = download(basis.x, p.dim * 6);
        const std::vector<Z> sx = download(basis.sx, p.dim * 6);
        const std::vector<Z> hx = download(basis.hx, p.dim * 6);
        for (int i = 0; i < p.dim * 4; ++i) { EXPECT_EQ(repaired[i], x[i]); }
        for (int j = 4; j < 6; ++j)
        {
            for (int row = 0; row < p.dim; ++row)
            {
                if (!p.physical(row)) { continue; }
                const int i = p.global_row(row);
                EXPECT_LT(std::abs(sx[row + j * p.dim] - (1 + .05 * i) * repaired[row + j * p.dim]), 1e-12);
                EXPECT_LT(std::abs(hx[row + j * p.dim] - (i + 1.0) * sx[row + j * p.dim]), 1e-12);
            }
            for (int i = 0; i < 6; ++i)
            {
                Z gram(0);
                for (int row = 0; row < p.dim; ++row) { gram += std::conj(repaired[row + i * p.dim]) * sx[row + j * p.dim]; }
                sum(&gram, 1);
                const Z expected = i == j ? Z(1) : Z(0);
                EXPECT_LT(std::abs(gram - expected), 1e-12);
            }
        }
        std::fill_n(x.data(), p.dim, Z(0));
        upload(basis.x, x);
        EXPECT_LT(solver.repair_basis(p, 3), 0);
    }

};
} // namespace hsolver

namespace
{
TEST_P(Lobpcg, StandardEigenpairs)
{
    Options o;
    o.overlap = false;
    check_case<Z>(o, GetParam());
    check_case<std::complex<float>>(o, GetParam());
}
TEST_P(Lobpcg, GeneralizedEigenpairs)
{
    check_case<Z>(Options{}, GetParam());
    check_case<std::complex<float>>(Options{}, GetParam());
}
TEST_P(Lobpcg, ExtraBands)
{
    Options o;
    o.extra = 5;
    check_case<Z>(o, GetParam());
    check_case<std::complex<float>>(o, GetParam());
}
TEST_P(Lobpcg, SpinorPadding)
{
    Options o;
    o.npol = 2;
    check_case<Z>(o, GetParam());
    check_case<std::complex<float>>(o, GetParam());
}
TEST_P(Lobpcg, DegenerateSpectrum)
{
    Options o;
    o.extra = 4;
    Problem<Z> p(o, GetParam(), comm_info());
    p.diagonal_problem();
    for (int i = 3; i < 8; ++i)
    {
        p.h[i + i * o.n] = 4.0 * p.s[i + i * o.n];
        p.reference_values[i] = 4;
    }
    ASSERT_TRUE(p.solve(1e-12, 250, {}));
    p.check(true, 1e-8, 1e-8);
}
TEST_P(Lobpcg, PartiallyExactInitialVectors)
{
    Problem<Z> p(Options{}, GetParam(), comm_info());
    p.exact_initial(2);
    ASSERT_TRUE(p.solve(1e-12, 250, {}));
    p.check(true, 1e-8, 1e-8);
}
TEST_P(Lobpcg, EigenvectorsMatchEigenvalues)
{
    check_eigenpair_identity<Z>(GetParam());
    check_eigenpair_identity<std::complex<float>>(GetParam());
}
TEST_P(Lobpcg, ExactInitialVectors)
{
    check_exact_initial<Z>(GetParam());
    check_exact_initial<std::complex<float>>(GetParam());
}
#ifdef __MPI
TEST_P(Lobpcg, DistributedRitz)
{
    Options o;
    o.backend = 2;
    check_case<Z>(o, GetParam());
    check_case<std::complex<float>>(o, GetParam());
}
#endif
TEST_P(Lobpcg, BudgetAndReuse)
{
    Options o;
    Problem<Z> p(o, GetParam(), comm_info());
    EXPECT_FALSE(p.solve(1e-13, 1, {}));
    EXPECT_EQ(p.solver->iterations(), 1);
    p.check(false, 1e-8, 1e-8);
    EXPECT_TRUE(p.solve(1e-13, 250, {}));
    p.check(true, 1e-8, 1e-8);
}
TEST(LobpcgInitialization, ExactNonorthogonalInitialization)
{
    Options o;
    o.bands = 4;
    Problem<Z> p(o, Locking::Hard, comm_info());
    p.exact_initial(o.bands);
    const std::vector<Z> original = p.input;
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
    EXPECT_TRUE(p.solve(1e-10, 250, {}));
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
    ASSERT_TRUE(implicit.solve(1e-12, 250, {}));
    EXPECT_EQ(implicit.s_calls, 0);
    o.overlap = true;
    Problem<Z> explicit_s(o, GetParam(), comm_info());
    explicit_s.s = implicit.s;
    explicit_s.reference_values = implicit.reference_values;
    ASSERT_TRUE(explicit_s.solve(1e-12, 250, {}));
    EXPECT_GT(explicit_s.s_calls, 0);
    explicit_s.check(true, 1e-8, 1e-8);
    for (int i = 0; i < o.bands; ++i)
    {
        EXPECT_NEAR(explicit_s.eigenvalues[i], implicit.eigenvalues[i], 1e-9);
    }
}
TEST_P(Lobpcg, BandThresholdOverridesScalar)
{
    Options o;
    o.extra = 4;
    Problem<Z> p(o, GetParam(), comm_info());
    const std::vector<double> thresholds(o.bands, 1e-13);
    ASSERT_TRUE(p.solve(1e-3, 250, thresholds));
    p.check(true, 1e-8, 1e-8);
}
TEST_P(Lobpcg, PrefixReactivationPolicy) { hsolver::LobpcgProjectionTest::prefix(GetParam()); }
TEST(LobpcgSubspace, DirectionImages) { hsolver::LobpcgProjectionTest::direction(0); }
TEST(LobpcgSubspace, DirectionRemovesLockedVector) { hsolver::LobpcgProjectionTest::direction(1); }
TEST(LobpcgSubspace, RetainedProjection) { hsolver::LobpcgProjectionTest::retained_projection(0); }
TEST(LobpcgSubspace, RetainedProjectionWithZeroPrefix) { hsolver::LobpcgProjectionTest::retained_projection(2); }
TEST(LobpcgSubspace, RefreshFailedMetric) { hsolver::LobpcgProjectionTest::recovery(); }
TEST(LobpcgSubspace, RepairDependentIncrement) { hsolver::LobpcgProjectionTest::rank_repair(); }
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
    testing::InitGoogleTest(&argc, argv);
    testing::TestEventListeners& listeners = testing::UnitTest::GetInstance()->listeners();
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
