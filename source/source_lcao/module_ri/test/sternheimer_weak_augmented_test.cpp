#include "source_lcao/module_ri/sternheimer_weak_augmented.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
using Blocks = ModuleRI::SternheimerWeakAugmented;
using Worker = Blocks::Worker;
using Complex = Blocks::Complex;
using Vector = Blocks::Vector;
using Matrix = std::vector<Vector>;

Complex dot(const Vector& a, const Vector& b)
{
    Complex value = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        value += std::conj(a[i]) * b[i];
    }
    return value;
}

Vector multiply(const Matrix& a, const Vector& x)
{
    Vector y(a.size(), 0.0);
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        for (std::size_t j = 0; j < x.size(); ++j)
        {
            y[i] += a[i][j] * x[j];
        }
    }
    return y;
}

// Independent small dense reference, deliberately not the helper's LAPACK path.
Vector direct_solve(Matrix a, Vector b)
{
    for (std::size_t k = 0; k < b.size(); ++k)
    {
        std::size_t pivot = k;
        for (std::size_t i = k + 1; i < b.size(); ++i)
        {
            if (std::abs(a[i][k]) > std::abs(a[pivot][k]))
            {
                pivot = i;
            }
        }
        if (std::abs(a[pivot][k]) < 1.0e-14)
        {
            throw std::runtime_error("singular test reference");
        }
        std::swap(a[k], a[pivot]);
        std::swap(b[k], b[pivot]);
        for (std::size_t i = k + 1; i < b.size(); ++i)
        {
            const Complex factor = a[i][k] / a[k][k];
            for (std::size_t j = k + 1; j < b.size(); ++j)
            {
                a[i][j] -= factor * a[k][j];
            }
            b[i] -= factor * b[k];
        }
    }
    for (std::size_t i = b.size(); i-- > 0;)
    {
        for (std::size_t j = i + 1; j < b.size(); ++j)
        {
            b[i] -= a[i][j] * b[j];
        }
        b[i] /= a[i][i];
    }
    return b;
}

void expect_vector(const Vector& actual, const Vector& expected, const double tolerance = 2.0e-11)
{
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
        EXPECT_NEAR(std::abs(actual[i] - expected[i]), 0.0, tolerance) << "entry " << i;
    }
}

Vector probe(const std::size_t n, const double offset)
{
    Vector v(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        v[i] = Complex(std::sin(0.8 * i + offset), std::cos(0.4 * i - offset));
    }
    return v;
}

struct Model
{
    static constexpr int nfine = 8;
    int nocc = 1;
    int nv = 2;
    int nc = 4;
    Matrix fine_h;
    Matrix e; // Columns in the independent fine Euclidean space.
    Matrix t; // Columns [F,W].
    Matrix h;
    Matrix s;
    Matrix ce;
    Blocks::Data data;

    Model()
    {
        fine_h.assign(nfine, Vector(nfine, 0.0));
        for (int i = 0; i < nfine; ++i)
        {
            fine_h[i][i] = 2.0 + 0.7 * i;
            for (int j = 0; j < i; ++j)
            {
                fine_h[i][j] = Complex(0.09 * std::cos(i + 2 * j), 0.07 * std::sin(2 * i + j));
                fine_h[j][i] = std::conj(fine_h[i][j]);
            }
        }
        for (int col = 0; col < nc; ++col)
        {
            Vector v(nfine);
            for (int i = 0; i < nfine; ++i)
            {
                v[i] = Complex(0.11 * std::sin((i + 1) * (col + 2)),
                               0.13 * std::cos((i + 2) * (col + 1)));
            }
            v[3 + col] += 1.0;
            for (const auto& previous : e)
            {
                const Complex overlap = dot(previous, v);
                for (int i = 0; i < nfine; ++i)
                {
                    v[i] -= overlap * previous[i];
                }
            }
            const double norm = std::sqrt(dot(v, v).real());
            for (auto& value : v)
            {
                value /= norm;
            }
            e.push_back(v);
        }
        const int nu = nocc + nv;
        data.nocc = nocc;
        data.nvirtual = nv;
        data.ncoarse = nc;
        data.hu.resize(nu * nu);
        data.l.resize(nu * nc);
        data.k.resize(nu * nc);
        for (int j = 0; j < nu; ++j)
        {
            for (int i = 0; i < nu; ++i)
            {
                data.hu[i + nu * j] = fine_h[i][j];
            }
        }
        ce.assign(nc, Vector(nc));
        for (int j = 0; j < nc; ++j)
        {
            const Vector he = multiply(fine_h, e[j]);
            for (int i = 0; i < nu; ++i)
            {
                data.l[i + nu * j] = e[j][i];
                data.k[i + nu * j] = he[i];
            }
            for (int i = 0; i < nc; ++i)
            {
                ce[i][j] = dot(e[i], he);
            }
        }
        for (int i = 0; i < nv; ++i)
        {
            Vector f(nfine, 0.0);
            f[nocc + i] = 1.0;
            t.push_back(f);
        }
        for (const auto& column : e)
        {
            Vector w = column;
            std::fill(w.begin(), w.begin() + nu, Complex(0.0));
            t.push_back(w);
        }
        h.assign(t.size(), Vector(t.size()));
        s = h;
        for (std::size_t j = 0; j < t.size(); ++j)
        {
            const Vector ht = multiply(fine_h, t[j]);
            for (std::size_t i = 0; i < t.size(); ++i)
            {
                h[i][j] = dot(t[i], ht);
                s[i][j] = dot(t[i], t[j]);
            }
        }
    }

    Blocks::Apply callback() const
    {
        const Matrix matrix = ce;
        return [matrix](const Vector& x, Vector& y) { y = multiply(matrix, x); };
    }

    Matrix shifted(const double eps, const double omega) const
    {
        Matrix result = h;
        for (std::size_t i = 0; i < h.size(); ++i)
        {
            for (std::size_t j = 0; j < h.size(); ++j)
            {
                result[i][j] -= Complex(eps, -omega) * s[i][j];
            }
        }
        return result;
    }
};

Blocks::Vertices vertices(const Model& model)
{
    Blocks::Vertices result;
    result.f = probe(model.nv, 0.3);
    result.w = probe(model.nc, -0.4);
    return result;
}

Vector joined(const Blocks::Vertices& g)
{
    Vector result = g.f;
    result.insert(result.end(), g.w.begin(), g.w.end());
    return result;
}

} // namespace

TEST(SternheimerWeakAugmented, MetricAndHamiltonianMatchIndependentFineSpace)
{
    const Model model;
    const auto blocks = std::make_shared<const Blocks>(model.data);
    Worker worker(blocks, model.callback(), 0.6, 0.2);
    const Vector v = probe(model.nv + model.nc, 0.8);
    Vector output;
    blocks->apply_metric(v, output);
    expect_vector(output, multiply(model.s, v));
    EXPECT_GT(std::abs(model.s[model.nv][model.nv] - 1.0), 0.01);
    EXPECT_GT(std::abs(model.s[model.nv][model.nv + 1]), 1.0e-4);
    EXPECT_GT(std::abs(model.h[0][model.nv].imag()), 1.0e-4);
    worker.apply_hamiltonian(v, output);
    expect_vector(output, multiply(model.h, v));
    worker.apply_shifted(v, output);
    expect_vector(output, multiply(model.shifted(0.6, 0.2), v));
    output = v;
    worker.apply_shifted(output, output);
    expect_vector(output, multiply(model.shifted(0.6, 0.2), v));
    output = v;
    blocks->apply_metric(output, output);
    expect_vector(output, multiply(model.s, v));
}

TEST(SternheimerWeakAugmented, HermitianAdjointsAndFrequencySignArePreserved)
{
    const Model model;
    const auto blocks = std::make_shared<const Blocks>(model.data);
    const Vector a = probe(model.h.size(), 0.2);
    const Vector b = probe(model.h.size(), -0.6);
    Worker positive(blocks, model.callback(), 0.6, 0.7);
    Worker negative(blocks, model.callback(), 0.6, -0.7);
    Vector ha, hb;
    positive.apply_hamiltonian(a, ha);
    positive.apply_hamiltonian(b, hb);
    EXPECT_NEAR(std::abs(dot(a, hb) - dot(ha, b)), 0.0, 1.0e-11);
    positive.apply_shifted(b, hb);
    negative.apply_shifted(a, ha);
    EXPECT_NEAR(std::abs(dot(a, hb) - dot(ha, b)), 0.0, 1.0e-11);
}

TEST(SternheimerWeakAugmented, SchurRhsReconstructionAndSolveMatchDenseAtMultipleFrequencies)
{
    const Model model;
    const auto blocks = std::make_shared<const Blocks>(model.data);
    const auto g = vertices(model);
    ModuleRI::SternheimerRPA::SolverOptions options;
    options.max_iter = 20;
    options.residual_tol = 1.0e-11;
    for (const double omega : {0.0, 0.03, 0.7, 3.0})
    {
        SCOPED_TRACE(omega);
        Worker worker(blocks, model.callback(), 0.6, omega);
        const Matrix shifted = model.shifted(0.6, omega);
        Matrix av(model.nv, Vector(model.nv));
        for (int i = 0; i < model.nv; ++i)
        {
            std::copy_n(shifted[i].begin(), model.nv, av[i].begin());
        }
        const Vector x = probe(model.nc, 0.1);
        Vector bx(model.nv, 0.0);
        Vector expected(model.nc, 0.0);
        for (int i = 0; i < model.nv; ++i)
        {
            for (int j = 0; j < model.nc; ++j)
            {
                bx[i] += shifted[i][model.nv + j] * x[j];
            }
        }
        const Vector inverse_bx = direct_solve(av, bx);
        const Vector inverse_g = direct_solve(av, g.f);
        Vector expected_rhs = g.w;
        for (int i = 0; i < model.nc; ++i)
        {
            expected_rhs[i] = -g.w[i];
            for (int j = 0; j < model.nc; ++j)
            {
                expected[i] += shifted[model.nv + i][model.nv + j] * x[j];
            }
            for (int j = 0; j < model.nv; ++j)
            {
                expected[i] -= shifted[model.nv + i][j] * inverse_bx[j];
                expected_rhs[i] += shifted[model.nv + i][j] * inverse_g[j];
            }
        }
        Vector actual;
        worker.apply_schur(x, actual);
        expect_vector(actual, expected);
        expect_vector(worker.schur_rhs(g), expected_rhs);
        Vector rhs = joined(g);
        for (auto& value : rhs)
        {
            value = -value;
        }
        const Vector direct = direct_solve(shifted, rhs);
        const Vector direct_x(direct.begin() + model.nv, direct.end());
        expect_vector(worker.reconstruct(g, direct_x), direct);
        const auto result = worker.solve(g, options, model.nc);
        ASSERT_TRUE(result.converged);
        EXPECT_TRUE(result.schur.converged);
        EXPECT_LE(result.relative_residual, options.residual_tol);
        expect_vector(result.coefficients, direct, 2.0e-10);
        const Complex response = dot(joined(g), result.coefficients);
        EXPECT_LT(response.real(), 0.0);
        if (omega > 0.0)
        {
            EXPECT_GT(response.imag(), 0.0); // -g*(H-eps*S+i*omega*S)^(-1)g.
        }
        else
        {
            EXPECT_NEAR(response.imag(), 0.0, 1.0e-11);
        }
    }
}

TEST(SternheimerWeakAugmented, ProjectsBothVerticesAndReconstructsFineResponseWithoutWeights)
{
    const Model model;
    const Blocks blocks(model.data);
    const Vector perturbation = probe(Model::nfine, 0.6);
    const Vector gu(perturbation.begin(), perturbation.begin() + model.nocc + model.nv);
    Vector ge;
    for (const auto& e : model.e)
    {
        ge.push_back(dot(e, perturbation));
    }
    Vector expected;
    for (const auto& t : model.t)
    {
        expected.push_back(dot(t, perturbation));
    }
    expect_vector(joined(blocks.project_vertices(gu, ge)), expected);
    const Vector coefficients = probe(model.t.size(), -0.2);
    const auto expansion = blocks.expand_coordinates(coefficients);
    Vector fine(Model::nfine, 0.0);
    Vector direct(Model::nfine, 0.0);
    std::copy(expansion.u.begin(), expansion.u.end(), fine.begin());
    for (int i = 0; i < Model::nfine; ++i)
    {
        for (int j = 0; j < model.nc; ++j)
        {
            fine[i] += model.e[j][i] * expansion.e[j];
        }
        for (std::size_t j = 0; j < model.t.size(); ++j)
        {
            direct[i] += model.t[j][i] * coefficients[j];
        }
    }
    expect_vector(fine, direct);
    EXPECT_NEAR(std::abs(fine[0]), 0.0, 1.0e-12);
}

TEST(SternheimerWeakAugmented, SingularIndefiniteAndNearlySingularMetricAreDiagnosed)
{
    Blocks::Data data;
    data.nvirtual = 1;
    data.ncoarse = 2;
    data.hu = {2.0};
    data.k = {0.3, Complex(0.0, 0.2)};
    // Neither column has norm one; their coherent combination is null in W.
    data.l = {0.6, Complex(0.0, 0.8)};
    try
    {
        const Blocks blocks(data);
        FAIL() << "singular metric accepted";
    }
    catch (const std::domain_error& error)
    {
        EXPECT_NE(std::string(error.what()).find("metric"), std::string::npos);
        EXPECT_NE(std::string(error.what()).find("singular"), std::string::npos);
    }
    data.l[0] = 0.7;
    EXPECT_THROW(Blocks{data}, std::domain_error);
    data.l = {std::sqrt(1.0 - 1.0e-13), 0.0};
    EXPECT_THROW(Blocks{data}, std::domain_error);
    data.l = {0.6, 0.0};
    const Blocks regular(data);
    EXPECT_NEAR(regular.minimum_metric_eigenvalue(), 0.64, 1.0e-13);
}

TEST(SternheimerWeakAugmented, RejectsSingularVirtualResolventWithoutEnergyShift)
{
    Blocks::Data data;
    data.nvirtual = 2;
    data.hu = {1.0, Complex(0.0, -1.0), Complex(0.0, 1.0), 1.0};
    const auto blocks = std::make_shared<const Blocks>(data);
    EXPECT_THROW((Worker(blocks, {}, 0.0, 0.0)), std::domain_error);
    Worker nonsingular(blocks, {}, 0.0, 0.2);
    const Blocks::Vertices g{{1.0, Complex(0.3, 0.2)}, {}};
    const auto result = nonsingular.solve(g, ModuleRI::SternheimerRPA::SolverOptions{}, 2);
    ASSERT_TRUE(result.converged);
    expect_vector(result.coefficients,
                  direct_solve({{Complex(1.0, 0.2), Complex(0.0, 1.0)},
                                {Complex(0.0, -1.0), Complex(1.0, 0.2)}},
                               {-g.f[0], -g.f[1]}));
}

TEST(SternheimerWeakAugmented, PureComplementAndZeroRightHandSideAreSupported)
{
    Blocks::Data data;
    data.ncoarse = 2;
    const auto blocks = std::make_shared<const Blocks>(data);
    const Matrix ce{{2.0, Complex(0.2, 0.1)}, {Complex(0.2, -0.1), 3.0}};
    Worker worker(blocks, [ce](const Vector& x, Vector& y) { y = multiply(ce, x); }, 0.5, 0.1);
    const auto zero = worker.solve({{}, {0.0, 0.0}}, ModuleRI::SternheimerRPA::SolverOptions{}, 2);
    EXPECT_TRUE(zero.converged);
    expect_vector(zero.coefficients, {0.0, 0.0});
    EXPECT_EQ(zero.absolute_residual, 0.0);
    const auto result = worker.solve({{}, {1.0, Complex(0.1, 0.3)}},
                                     ModuleRI::SternheimerRPA::SolverOptions{}, 2);
    ASSERT_TRUE(result.converged);
    expect_vector(result.coefficients,
                  direct_solve({{Complex(1.5, 0.1), Complex(0.2, 0.1)},
                                {Complex(0.2, -0.1), Complex(2.5, 0.1)}},
                               {-1.0, Complex(-0.1, -0.3)}));
}

TEST(SternheimerWeakAugmented, NumericallyUnresolvedVirtualResolventIsDiagnosed)
{
    Blocks::Data data;
    data.nvirtual = 2;
    data.hu = {1.0, 0.0, 0.0, 1.0e-20};
    const auto blocks = std::make_shared<const Blocks>(data);
    EXPECT_THROW((Worker(blocks, {}, 0.0, 0.0)), std::domain_error);
    // A genuinely small but resolved energy is not replaced by a shifted one.
    data.hu[3] = 1.0e-8;
    const auto resolved = std::make_shared<const Blocks>(data);
    Worker worker(resolved, {}, 0.0, 0.0);
    expect_vector(worker.reconstruct({{0.0, 1.0e-8}, {}}, {}), {0.0, -1.0});
}

TEST(SternheimerWeakAugmented, NonconvergenceReturnsActualFullBlockResidual)
{
    const Model model;
    const auto blocks = std::make_shared<const Blocks>(model.data);
    Worker worker(blocks, model.callback(), 0.6, 0.3);
    ModuleRI::SternheimerRPA::SolverOptions options;
    options.max_iter = 1;
    options.residual_tol = 1.0e-14;
    const auto g = vertices(model);
    const auto result = worker.solve(g, options, 1);
    EXPECT_FALSE(result.converged);
    EXPECT_FALSE(result.schur.converged);
    Vector residual = multiply(model.shifted(0.6, 0.3), result.coefficients);
    const Vector source = joined(g);
    for (std::size_t i = 0; i < residual.size(); ++i)
    {
        residual[i] += source[i];
    }
    EXPECT_NEAR(result.absolute_residual, std::sqrt(dot(residual, residual).real()), 1.0e-12);
    EXPECT_NEAR(result.relative_residual,
                std::sqrt(dot(residual, residual).real() / dot(source, source).real()), 1.0e-12);
}

TEST(SternheimerWeakAugmented, RejectsMalformedDataVectorsCallbacksAndNonfiniteValues)
{
    const Model model;
    auto bad = model.data;
    bad.l.pop_back();
    EXPECT_THROW(Blocks{bad}, std::invalid_argument);
    bad = model.data;
    bad.k[0] = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(Blocks{bad}, std::invalid_argument);
    bad = model.data;
    bad.hu[1] += Complex(0.0, 0.5);
    EXPECT_THROW(Blocks{bad}, std::invalid_argument);
    bad = model.data;
    bad.nocc = -1;
    EXPECT_THROW(Blocks{bad}, std::invalid_argument);
    EXPECT_THROW(Blocks{Blocks::Data{}}, std::invalid_argument);
    const auto blocks = std::make_shared<const Blocks>(model.data);
    EXPECT_THROW((Worker(nullptr, model.callback(), 0.6, 0.1)), std::invalid_argument);
    EXPECT_THROW((Worker(blocks, {}, 0.6, 0.1)), std::invalid_argument);
    EXPECT_THROW((Worker(blocks, model.callback(), std::numeric_limits<double>::infinity(), 0.1)),
                 std::invalid_argument);
    Worker worker(blocks, model.callback(), 0.6, 0.1);
    Vector output;
    EXPECT_THROW(worker.apply_shifted({}, output), std::invalid_argument);
    EXPECT_THROW(worker.schur_rhs({{}, {}}), std::invalid_argument);
    EXPECT_THROW(blocks->project_vertices({}, {}), std::invalid_argument);
    Vector bad_vector = probe(model.nc, 0.2);
    bad_vector[0] = Complex(0.0, std::numeric_limits<double>::infinity());
    EXPECT_THROW(worker.apply_schur(bad_vector, output), std::invalid_argument);
    Worker wrong_size(blocks, [](const Vector&, Vector& y) { y.assign(1, 0.0); }, 0.6, 0.1);
    EXPECT_THROW(wrong_size.apply_schur(probe(model.nc, 0.2), output), std::invalid_argument);
    Worker nonfinite(blocks, [](const Vector& x, Vector& y) {
        y.assign(x.size(), std::numeric_limits<double>::quiet_NaN());
    }, 0.6, 0.1);
    EXPECT_THROW(nonfinite.apply_schur(probe(model.nc, 0.2), output), std::invalid_argument);
    ModuleRI::SternheimerRPA::SolverOptions invalid;
    invalid.residual_tol = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(worker.solve(vertices(model), invalid, 3), std::invalid_argument);
    EXPECT_THROW(worker.solve(vertices(model), ModuleRI::SternheimerRPA::SolverOptions{}, 0),
                 std::invalid_argument);
}

TEST(SternheimerWeakAugmented, EnforcesPayloadBoundsAndOwnsImmutableInputCopies)
{
    Model model;
    EXPECT_THROW((Blocks(model.data, 1)), std::length_error);
    const auto blocks = std::make_shared<const Blocks>(model.data);
    EXPECT_THROW((Worker(blocks, model.callback(), 0.6, 0.1, 1)), std::length_error);
    Worker first(blocks, model.callback(), 0.6, 0.1);
    Worker second(blocks, model.callback(), 0.6, 0.8);
    const Vector v = probe(model.h.size(), 0.4);
    model.data.hu[0] += 10.0;
    model.data.l[0] += 0.3;
    Vector a, b;
    first.apply_shifted(v, a);
    second.apply_shifted(v, b);
    expect_vector(a, multiply(model.shifted(0.6, 0.1), v));
    expect_vector(b, multiply(model.shifted(0.6, 0.8), v));
    first.apply_shifted(v, b);
    expect_vector(a, b);
    const auto required = first.workspace_bytes_required(3);
    Worker bounded(blocks, model.callback(), 0.6, 0.1, required);
    EXPECT_THROW(bounded.solve(vertices(model), ModuleRI::SternheimerRPA::SolverOptions{}, 1000),
                 std::length_error);
}

TEST(SternheimerWeakAugmented, LargeCoarseSpaceDoesNotRequireDenseCoarseMatrices)
{
    Blocks::Data data;
    data.nvirtual = 1;
    data.ncoarse = 10000;
    data.hu = {2.0};
    data.l.assign(data.ncoarse, 0.001);
    data.k.assign(data.ncoarse, Complex(0.002, 0.001));
    const auto blocks = std::make_shared<const Blocks>(data, 1024 * 1024);
    EXPECT_LT(blocks->storage_bytes(), 1024U * 1024);
    EXPECT_NEAR(blocks->minimum_metric_eigenvalue(), 0.99, 1.0e-12);
    Worker worker(blocks, [](const Vector& x, Vector& y) {
        y = x;
        for (auto& value : y)
        {
            value *= 3.0;
        }
    }, 0.5, 0.2, 8 * 1024 * 1024);
    Vector result;
    worker.apply_schur(Vector(data.ncoarse, 1.0), result);
    ASSERT_EQ(result.size(), static_cast<std::size_t>(data.ncoarse));
    const Complex expected = 2.98 - Complex(0.5, -0.2) * 0.99 - 0.01 / Complex(1.5, 0.2);
    EXPECT_NEAR(std::abs(result[0] - expected), 0.0, 1.0e-10);
    EXPECT_NEAR(std::abs(result.back() - expected), 0.0, 1.0e-10);
}

TEST(SternheimerWeakAugmented, ResolvesSmallPositiveMetricWithoutDroppingModes)
{
    const int n = 32;
    Blocks::Data data;
    data.nocc = n;
    data.nvirtual = 0;
    data.ncoarse = n;
    data.hu.assign(n * n, 0.0);
    data.k.assign(n * n, 0.0);
    data.l.assign(n * n, 0.0);
    Vector metric(n), rhs(n), expected(n);
    for (int i = 0; i < n; ++i)
    {
        const double l = std::sqrt(1.0 - std::pow(10.0, -8.0 * i / (n - 1)));
        data.l[i + n * i] = l;
        metric[i] = 1.0 - l * l;
        rhs[i] = std::sqrt(metric[i].real());
        expected[i] = -rhs[i] / (Complex(2.0, 0.2) * metric[i]);
    }
    auto blocks = std::make_shared<Blocks>(std::move(data));
    Worker worker(blocks, [metric](const Vector& x, Vector& out) {
        out.resize(x.size());
        for (std::size_t i = 0; i < x.size(); ++i) out[i] = 2.0 * metric[i] * x[i];
    }, 0.0, 0.2);
    ModuleRI::SternheimerRPA::SolverOptions options;
    options.max_iter = 40;
    options.residual_tol = 1e-8;
    const auto result = worker.solve({{}, rhs}, options, 8);
    ASSERT_TRUE(result.converged) << result.relative_residual;
    EXPECT_LE(result.relative_residual, options.residual_tol);
    EXPECT_LE(result.schur.iterations, 5);
    expect_vector(result.coefficients, expected, 1e-4);
    EXPECT_GT(blocks->minimum_metric_eigenvalue(), 9e-9);
}

TEST(SternheimerWeakAugmented, NormalizedStoppingPreservesOriginalResidualTolerance)
{
    Blocks::Data data;
    data.nocc = 1;
    data.ncoarse = 2;
    data.hu = {0.0};
    const double l = std::sqrt(1.0 - 1e-8);
    const double s = 1.0 - l * l;
    data.l = {l, 0.0};
    data.k = {0.0, 0.0};
    auto blocks = std::make_shared<Blocks>(data);
    Worker worker(blocks, [s](const Vector& x, Vector& out) {
        out = {s * x[0], 2.0 * x[1]};
    }, 0.0, 0.0);
    ModuleRI::SternheimerRPA::SolverOptions options;
    options.max_iter = 10;
    options.residual_tol = 1e-8;
    const auto result = worker.solve({{}, {-1.0, -1e-6}}, options, 2);
    ASSERT_TRUE(result.converged) << result.relative_residual;
    EXPECT_LE(result.relative_residual, options.residual_tol);
    EXPECT_NEAR(std::abs(s * result.coefficients[0] - Complex(1.0)), 0.0, 1e-8);
    EXPECT_NEAR(std::abs(result.coefficients[1] - Complex(5e-7)), 0.0, 1e-9);
}

TEST(SternheimerWeakAugmented, InverseMetricRootNormalizesComplexOverlapAndPreservesAdjoint)
{
    Model model;
    Blocks blocks(model.data);
    const Vector x = probe(model.nc, 0.6), y = probe(model.nc, 1.1);
    Vector tx, ty;
    blocks.apply_complement_inverse_sqrt(x, tx);
    blocks.apply_complement_inverse_sqrt(y, ty);
    Vector augmented(model.nv, 0.0), metric;
    augmented.insert(augmented.end(), tx.begin(), tx.end());
    blocks.apply_metric(augmented, metric);
    Vector lower(metric.begin() + model.nv, metric.end()), normalized;
    blocks.apply_complement_inverse_sqrt(lower, normalized);
    expect_vector(normalized, x);
    EXPECT_NEAR(std::abs(dot(x, ty) - dot(tx, y)), 0.0, 2e-11);
    EXPECT_GT(dot(x, tx).real(), 0.0);
    Vector alias = x;
    blocks.apply_complement_inverse_sqrt(alias, alias);
    expect_vector(alias, tx);
    EXPECT_THROW(blocks.apply_complement_inverse_sqrt(Vector(1), normalized), std::invalid_argument);
    Vector invalid = x;
    invalid[0] = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(blocks.apply_complement_inverse_sqrt(invalid, normalized), std::invalid_argument);
}

TEST(SternheimerWeakAugmented, InverseMetricRootHandlesIdentityAndEmptyComplement)
{
    Blocks::Data data;
    data.ncoarse = 3;
    Blocks identity(data);
    Vector x = probe(3, 0.3), result;
    identity.apply_complement_inverse_sqrt(x, result);
    expect_vector(result, x);
    data.ncoarse = 0;
    data.nvirtual = 1;
    data.hu = {2.0};
    Blocks empty(data);
    empty.apply_complement_inverse_sqrt({}, result);
    EXPECT_TRUE(result.empty());
}

TEST(SternheimerWeakAugmented, PositiveVirtualBlockDoesNotImplyPositiveAugmentedHamiltonian)
{
    Blocks::Data data;
    data.nvirtual = 1;
    data.ncoarse = 1;
    data.hu = {2.0};
    data.l = {0.0};
    data.k = {Complex(3.0, 1.0)};
    const auto blocks = std::make_shared<const Blocks>(data);
    Worker worker(blocks, [](const Vector& x, Vector& y) { y = x; }, 0.0, 0.2);
    const Vector v{1.0, -1.0};
    Vector hv;
    worker.apply_hamiltonian(v, hv);
    EXPECT_LT(dot(v, hv).real(), 0.0);
    const auto result = worker.solve({{1.0}, {0.4}}, ModuleRI::SternheimerRPA::SolverOptions{}, 1);
    ASSERT_TRUE(result.converged);
    expect_vector(result.coefficients,
                  direct_solve({{Complex(2.0, 0.2), Complex(3.0, 1.0)},
                                {Complex(3.0, -1.0), Complex(1.0, 0.2)}}, {-1.0, -0.4}));
}
