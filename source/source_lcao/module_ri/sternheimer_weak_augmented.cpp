#include "source_lcao/module_ri/sternheimer_weak_augmented.h"

#include "source_base/module_external/blas_connector.h"
#include "source_base/module_external/lapack_connector.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

// The existing connector declares ZGETRF/ZHEEV but not these LU operations.
extern "C" void zgetrs_(const char* trans, const int* n, const int* nrhs,
                        const std::complex<double>* a, const int* lda, const int* ipiv,
                        std::complex<double>* b, const int* ldb, int* info);
extern "C" void zgecon_(const char* norm, const int* n, const std::complex<double>* a,
                        const int* lda, const double* anorm, double* rcond,
                        std::complex<double>* work, double* rwork, int* info);

namespace ModuleRI
{
namespace
{
using Blocks = SternheimerWeakAugmented;
using Complex = Blocks::Complex;
using Vector = Blocks::Vector;

std::size_t add(const std::size_t a, const std::size_t b)
{
    if (b > std::numeric_limits<std::size_t>::max() - a)
    {
        throw std::length_error("weak augmented payload size overflow");
    }
    return a + b;
}

std::size_t product(const std::size_t a, const std::size_t b)
{
    if (a != 0 && b > std::numeric_limits<std::size_t>::max() / a)
    {
        throw std::length_error("weak augmented payload size overflow");
    }
    return a * b;
}

bool finite(const Complex value)
{
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}

void validate(const Vector& values, const std::size_t size, const char* name)
{
    if (values.size() != size)
    {
        throw std::invalid_argument(std::string("weak augmented ") + name + " dimension mismatch");
    }
    for (const auto value : values)
    {
        if (!finite(value))
        {
            throw std::invalid_argument(std::string("weak augmented ") + name + " contains nonfinite values");
        }
    }
}

void budget(const std::size_t required, const std::size_t limit)
{
    if (required > limit)
    {
        throw std::length_error("weak augmented numerical payload exceeds budget");
    }
}

// Packed column-major matrix-vector product; never forms a coarse square block.
Vector multiply(const Vector& matrix, const int rows, const int cols, const Vector& x,
                const bool adjoint = false)
{
    Vector y(adjoint ? cols : rows, 0.0);
    if (rows != 0 && cols != 0)
    {
        BlasConnector::gemv(adjoint ? 'C' : 'N', rows, cols, Complex(1.0),
                            matrix.data(), rows, x.data(), 1, Complex(0.0), y.data(), 1);
    }
    validate(y, y.size(), "matrix product");
    return y;
}

double norm(const Vector& x)
{
    double value = 0.0;
    for (const auto element : x)
    {
        value = std::hypot(value, std::abs(element));
    }
    if (!std::isfinite(value))
    {
        throw std::runtime_error("weak augmented norm overflow");
    }
    return value;
}

Complex euclidean_dot(const Vector& a, const Vector& b)
{
    validate(a, b.size(), "dot input");
    validate(b, a.size(), "dot input");
    Complex value = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        value += std::conj(a[i]) * b[i];
    }
    if (!finite(value))
    {
        throw std::runtime_error("weak augmented dot overflow");
    }
    return value;
}

} // namespace

SternheimerWeakAugmented::SternheimerWeakAugmented(Data data,
                                                 const std::size_t max_payload_bytes,
                                                 const double metric_tolerance,
                                                 const double hermitian_tolerance)
    : data_(std::move(data))
{
    if (data_.nocc < 0 || data_.nvirtual < 0 || data_.ncoarse < 0
        || data_.nocc > std::numeric_limits<int>::max() / 4 - data_.nvirtual
        || data_.ncoarse > std::numeric_limits<int>::max() - data_.nvirtual
        || (data_.nvirtual == 0 && data_.ncoarse == 0))
    {
        throw std::invalid_argument("weak augmented invalid dimensions");
    }
    if (!std::isfinite(metric_tolerance) || metric_tolerance <= 0.0 || metric_tolerance >= 1.0
        || !std::isfinite(hermitian_tolerance) || hermitian_tolerance <= 0.0
        || hermitian_tolerance >= 1.0)
    {
        throw std::invalid_argument("weak augmented invalid diagnostic tolerance");
    }
    nu_ = data_.nocc + data_.nvirtual;
    const std::size_t square = product(nu_, nu_);
    const std::size_t rectangular = product(nu_, data_.ncoarse);
    validate(data_.hu, square, "HU");
    validate(data_.l, rectangular, "L");
    validate(data_.k, rectangular, "K");
    // ZHEEV minimum work: max(1,2*nu-1) complex, max(1,3*nu-2) real.
    const std::size_t scratch = add(product(add(square, product(2, nu_)), sizeof(Complex)),
                                    product(product(4, nu_), sizeof(double)));
    // In addition to the eigensolver matrix, retain the small inverse-root
    // kernel and a temporary weighted eigenvector matrix.
    budget(add(storage_bytes(), add(scratch, product(product(2, square), sizeof(Complex)))), max_payload_bytes);
    double scale = 1.0;
    for (const auto value : data_.hu)
    {
        scale = std::max(scale, std::abs(value));
    }
    if (!std::isfinite(scale))
    {
        throw std::invalid_argument("weak augmented HU magnitude overflow");
    }
    for (int j = 0; j < nu_; ++j)
    {
        for (int i = 0; i < nu_; ++i)
        {
            const auto ij = i + static_cast<std::size_t>(nu_) * j;
            const auto ji = j + static_cast<std::size_t>(nu_) * i;
            if (std::abs(data_.hu[ij] / scale - std::conj(data_.hu[ji]) / scale) > hermitian_tolerance)
            {
                throw std::invalid_argument("weak augmented HU is not Hermitian");
            }
        }
    }
    if (nu_ == 0 || data_.ncoarse == 0)
    {
        return;
    }

    // LL* and L*L share their nonzero spectrum, so this detects coherent null
    // directions of W even when every individual W column has nonzero norm.
    Vector gram(square, 0.0);
    BlasConnector::gemm_cm('N', 'C', nu_, nu_, data_.ncoarse, Complex(1.0),
                           data_.l.data(), nu_, data_.l.data(), nu_, Complex(0.0), gram.data(), nu_);
    validate(gram, square, "LL* metric diagnostic");
    std::vector<double> eigenvalues(nu_);
    std::vector<double> rwork(std::max(1, 3 * nu_ - 2));
    const int lwork = std::max(1, 2 * nu_ - 1);
    Vector work(lwork);
    const char jobz = 'V';
    const char uplo = 'U';
    int info = 0;
    zheev_(&jobz, &uplo, &nu_, gram.data(), &nu_, eigenvalues.data(), work.data(), &lwork,
           rwork.data(), &info);
    if (info != 0 || !std::isfinite(eigenvalues.back()))
    {
        throw std::runtime_error("weak augmented metric eigensolver failed");
    }
    minimum_metric_eigenvalue_ = 1.0 - eigenvalues.back();
    if (minimum_metric_eigenvalue_ <= metric_tolerance)
    {
        std::ostringstream message;
        message.precision(17);
        message << "weak augmented complement metric is "
                << (minimum_metric_eigenvalue_ < -metric_tolerance ? "indefinite" : "singular or ill-conditioned")
                << "; minimum eigenvalue=" << minimum_metric_eigenvalue_
                << ", rejection threshold=" << metric_tolerance << "; no regularization applied";
        throw std::domain_error(message.str());
    }
    Vector scaled = gram;
    for (int j = 0; j < nu_; ++j)
    {
        if (!std::isfinite(eigenvalues[j]))
            throw std::runtime_error("weak augmented metric spectrum is nonfinite");
        const double root = std::sqrt(1.0 - eigenvalues[j]);
        // Stable even for lambda=0: ((1-lambda)^(-1/2)-1)/lambda.
        const double factor = 1.0 / (root * (1.0 + root));
        for (int i = 0; i < nu_; ++i) scaled[i + static_cast<std::size_t>(nu_) * j] *= factor;
    }
    inverse_sqrt_kernel_.resize(square);
    BlasConnector::gemm_cm('N', 'C', nu_, nu_, nu_, Complex(1.0),
                           scaled.data(), nu_, gram.data(), nu_, Complex(0.0),
                           inverse_sqrt_kernel_.data(), nu_);
    validate(inverse_sqrt_kernel_, square, "metric inverse square root kernel");
}

int SternheimerWeakAugmented::nocc() const { return data_.nocc; }
int SternheimerWeakAugmented::nvirtual() const { return data_.nvirtual; }
int SternheimerWeakAugmented::ncoarse() const { return data_.ncoarse; }

std::size_t SternheimerWeakAugmented::storage_bytes() const
{
    return product(add(inverse_sqrt_kernel_.capacity(),
                       add(data_.hu.capacity(), add(data_.l.capacity(), data_.k.capacity()))), sizeof(Complex));
}

double SternheimerWeakAugmented::minimum_metric_eigenvalue() const
{
    return minimum_metric_eigenvalue_;
}

void SternheimerWeakAugmented::validate_coefficients(const Vector& coefficients) const
{
    validate(coefficients, add(nvirtual(), ncoarse()), "[y;x]");
}

void SternheimerWeakAugmented::apply_complement_metric(const Vector& x, Vector& output) const
{
    const Vector lx = multiply(data_.l, nu_, ncoarse(), x);
    Vector result = multiply(data_.l, nu_, ncoarse(), lx, true);
    for (int i = 0; i < ncoarse(); ++i)
    {
        result[i] = x[i] - result[i];
    }
    validate(result, ncoarse(), "S_W x");
    output = std::move(result);
}

void SternheimerWeakAugmented::apply_metric(const Vector& coefficients, Vector& output) const
{
    validate_coefficients(coefficients);
    const Vector x(coefficients.begin() + nvirtual(), coefficients.end());
    Vector sx;
    apply_complement_metric(x, sx);
    Vector result(coefficients.begin(), coefficients.begin() + nvirtual());
    result.insert(result.end(), sx.begin(), sx.end());
    output = std::move(result);
}

void SternheimerWeakAugmented::apply_complement_inverse_sqrt(const Vector& x, Vector& output) const
{
    validate(x, ncoarse(), "inverse metric square root input");
    if (nu_ == 0 || ncoarse() == 0)
    {
        output = x;
        return;
    }
    const Vector lx = multiply(data_.l, nu_, ncoarse(), x);
    const Vector correction = multiply(inverse_sqrt_kernel_, nu_, nu_, lx);
    Vector result = multiply(data_.l, nu_, ncoarse(), correction, true);
    for (int i = 0; i < ncoarse(); ++i) result[i] += x[i];
    validate(result, ncoarse(), "inverse metric square root output");
    output = std::move(result);
}

void SternheimerWeakAugmented::apply_b(const Vector& x, Vector& output) const
{
    const Vector lx = multiply(data_.l, nu_, ncoarse(), x);
    const Vector hlx = multiply(data_.hu, nu_, nu_, lx);
    const Vector kx = multiply(data_.k, nu_, ncoarse(), x);
    Vector result(nvirtual());
    for (int i = 0; i < nvirtual(); ++i)
    {
        result[i] = kx[nocc() + i] - hlx[nocc() + i];
    }
    validate(result, nvirtual(), "B x");
    output = std::move(result);
}

void SternheimerWeakAugmented::apply_b_adjoint(const Vector& y, Vector& output) const
{
    Vector uy(nu_, 0.0);
    std::copy(y.begin(), y.end(), uy.begin() + nocc());
    const Vector huy = multiply(data_.hu, nu_, nu_, uy, true);
    Vector result = multiply(data_.k, nu_, ncoarse(), uy, true);
    const Vector lhuy = multiply(data_.l, nu_, ncoarse(), huy, true);
    for (int i = 0; i < ncoarse(); ++i)
    {
        result[i] -= lhuy[i];
    }
    validate(result, ncoarse(), "B* y");
    output = std::move(result);
}

SternheimerWeakAugmented::Vertices SternheimerWeakAugmented::project_vertices(const Vector& g_u,
                                                                            const Vector& g_e) const
{
    validate(g_u, nu_, "gU");
    validate(g_e, ncoarse(), "gE");
    Vertices result;
    result.f.assign(g_u.begin() + nocc(), g_u.end());
    result.w = multiply(data_.l, nu_, ncoarse(), g_u, true);
    for (int i = 0; i < ncoarse(); ++i)
    {
        result.w[i] = g_e[i] - result.w[i];
    }
    validate(result.w, ncoarse(), "gW");
    return result;
}

SternheimerWeakAugmented::Expansion SternheimerWeakAugmented::expand_coordinates(const Vector& coefficients) const
{
    validate_coefficients(coefficients);
    Expansion result;
    result.e.assign(coefficients.begin() + nvirtual(), coefficients.end());
    result.u = multiply(data_.l, nu_, ncoarse(), result.e);
    for (auto& value : result.u)
    {
        value = -value;
    }
    for (int i = 0; i < nvirtual(); ++i)
    {
        result.u[nocc() + i] += coefficients[i];
    }
    validate(result.u, nu_, "expanded U coordinates");
    return result;
}

SternheimerWeakAugmented::Worker::Worker(std::shared_ptr<const SternheimerWeakAugmented> blocks,
                                        Apply ce, const double eps_source, const double omega,
                                        const std::size_t max_workspace_bytes)
    : blocks_(std::move(blocks)), ce_(std::move(ce)), z_(eps_source, -omega),
      max_workspace_bytes_(max_workspace_bytes)
{
    if (!blocks_ || (blocks_->ncoarse() != 0 && !ce_) || !finite(z_))
    {
        throw std::invalid_argument("weak augmented worker requires blocks, finite shift and CE callback");
    }
    budget(workspace_bytes_required(0), max_workspace_bytes_);
    const int nv = blocks_->nvirtual();
    if (nv == 0)
    {
        return;
    }
    lu_.resize(product(nv, nv));
    pivots_.resize(nv);
    for (int j = 0; j < nv; ++j)
    {
        for (int i = 0; i < nv; ++i)
        {
            lu_[i + static_cast<std::size_t>(nv) * j]
                = blocks_->data_.hu[blocks_->nocc() + i
                                   + static_cast<std::size_t>(blocks_->nu_) * (blocks_->nocc() + j)]
                  - (i == j ? z_ : Complex(0.0));
        }
    }
    validate(lu_, product(nv, nv), "virtual resolvent");
    double one_norm = 0.0;
    for (int j = 0; j < nv; ++j)
    {
        double column_sum = 0.0;
        for (int i = 0; i < nv; ++i)
        {
            column_sum += std::abs(lu_[i + static_cast<std::size_t>(nv) * j]);
        }
        one_norm = std::max(one_norm, column_sum);
    }
    if (!std::isfinite(one_norm))
    {
        throw std::runtime_error("weak augmented virtual resolvent norm overflow");
    }
    int info = 0;
    zgetrf_(&nv, &nv, lu_.data(), &nv, pivots_.data(), &info);
    if (info > 0)
    {
        throw std::domain_error("weak augmented singular virtual resolvent A-zI; no energy shift applied");
    }
    if (info < 0)
    {
        throw std::runtime_error("weak augmented virtual LU factorization failed");
    }
    validate(lu_, product(nv, nv), "virtual LU");
    const char norm_type = '1';
    double rcond = 0.0;
    Vector work(2 * nv);
    std::vector<double> rwork(2 * nv);
    zgecon_(&norm_type, &nv, lu_.data(), &nv, &one_norm, &rcond, work.data(), rwork.data(), &info);
    if (info != 0 || !std::isfinite(rcond))
    {
        throw std::runtime_error("weak augmented virtual resolvent condition estimate failed");
    }
    if (rcond <= std::numeric_limits<double>::epsilon())
    {
        std::ostringstream message;
        message.precision(17);
        message << "weak augmented singular or numerically unresolved virtual resolvent A-zI; rcond="
                << rcond << "; no energy shift applied";
        throw std::domain_error(message.str());
    }
}

std::size_t SternheimerWeakAugmented::Worker::workspace_bytes_required(const int restart_dimension) const
{
    if (restart_dimension < 0)
    {
        throw std::invalid_argument("weak augmented negative restart dimension");
    }
    const std::size_t nv = blocks_->nvirtual();
    const std::size_t nc = blocks_->ncoarse();
    const std::size_t nu = blocks_->nu_;
    const std::size_t r = restart_dimension;
    // Generous vector allowance covers nested products, full residual, RHS,
    // reconstruction and callback output (not opaque callback internals).
    std::size_t count = add(product(nv, nv), product(32, add(nu, nc)));
    count = add(count, product(add(product(2, r), 1), nc));
    // Hessenberg plus the unmodified GMRES normal-equation solver and copies.
    count = add(count, add(product(4, product(add(r, 1), add(r, 1))), product(8, add(r, 1))));
    return add(product(count, sizeof(Complex)), product(nv, sizeof(int)));
}

Vector SternheimerWeakAugmented::Worker::solve_virtual(const Vector& rhs) const
{
    const int n = blocks_->nvirtual();
    validate(rhs, n, "virtual RHS");
    Vector result = rhs;
    if (n != 0)
    {
        const char trans = 'N';
        const int nrhs = 1;
        int info = 0;
        zgetrs_(&trans, &n, &nrhs, lu_.data(), &n, pivots_.data(), result.data(), &n, &info);
        if (info != 0)
        {
            throw std::runtime_error("weak augmented virtual LU solve failed");
        }
        validate(result, n, "virtual solution");
    }
    return result;
}

void SternheimerWeakAugmented::Worker::apply_complement(const Vector& x, Vector& output) const
{
    const auto& d = blocks_->data_;
    const int nu = blocks_->nu_;
    Vector result(d.ncoarse, 0.0);
    if (d.ncoarse != 0)
    {
        ce_(x, result);
        validate(result, d.ncoarse, "CE output");
    }
    const Vector lx = multiply(d.l, nu, d.ncoarse, x);
    const Vector kx = multiply(d.k, nu, d.ncoarse, x);
    const Vector hlx = multiply(d.hu, nu, nu, lx);
    const Vector lkx = multiply(d.l, nu, d.ncoarse, kx, true);
    const Vector klx = multiply(d.k, nu, d.ncoarse, lx, true);
    const Vector lhlx = multiply(d.l, nu, d.ncoarse, hlx, true);
    for (int i = 0; i < d.ncoarse; ++i)
    {
        result[i] += -lkx[i] - klx[i] + lhlx[i];
    }
    validate(result, d.ncoarse, "C x");
    output = std::move(result);
}

void SternheimerWeakAugmented::Worker::apply_hamiltonian(const Vector& coefficients, Vector& output) const
{
    blocks_->validate_coefficients(coefficients);
    const int nv = blocks_->nvirtual();
    const Vector y(coefficients.begin(), coefficients.begin() + nv);
    const Vector x(coefficients.begin() + nv, coefficients.end());
    Vector result;
    blocks_->apply_b(x, result);
    for (int j = 0; j < nv; ++j)
    {
        for (int i = 0; i < nv; ++i)
        {
            result[i] += blocks_->data_.hu[blocks_->nocc() + i
                                          + static_cast<std::size_t>(blocks_->nu_) * (blocks_->nocc() + j)] * y[j];
        }
    }
    Vector lower;
    apply_complement(x, lower);
    Vector by;
    blocks_->apply_b_adjoint(y, by);
    for (int i = 0; i < blocks_->ncoarse(); ++i)
    {
        lower[i] += by[i];
    }
    result.insert(result.end(), lower.begin(), lower.end());
    validate(result, coefficients.size(), "H augmented output");
    output = std::move(result);
}

void SternheimerWeakAugmented::Worker::apply_shifted(const Vector& coefficients, Vector& output) const
{
    Vector h, s;
    apply_hamiltonian(coefficients, h);
    blocks_->apply_metric(coefficients, s);
    for (std::size_t i = 0; i < h.size(); ++i)
    {
        h[i] -= z_ * s[i];
    }
    validate(h, coefficients.size(), "(H-zS) output");
    output = std::move(h);
}

void SternheimerWeakAugmented::Worker::apply_schur(const Vector& x, Vector& output) const
{
    validate(x, blocks_->ncoarse(), "Schur input");
    Vector c, s, bx, correction;
    apply_complement(x, c);
    blocks_->apply_complement_metric(x, s);
    blocks_->apply_b(x, bx);
    blocks_->apply_b_adjoint(solve_virtual(bx), correction);
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        c[i] -= z_ * s[i] + correction[i];
    }
    validate(c, x.size(), "Schur output");
    output = std::move(c);
}

void SternheimerWeakAugmented::Worker::validate_vertices(const Vertices& g) const
{
    validate(g.f, blocks_->nvirtual(), "gF");
    validate(g.w, blocks_->ncoarse(), "gW");
}

Vector SternheimerWeakAugmented::Worker::schur_rhs(const Vertices& g) const
{
    validate_vertices(g);
    Vector result;
    blocks_->apply_b_adjoint(solve_virtual(g.f), result);
    for (std::size_t i = 0; i < result.size(); ++i)
    {
        result[i] -= g.w[i];
    }
    validate(result, blocks_->ncoarse(), "Schur RHS");
    return result;
}

Vector SternheimerWeakAugmented::Worker::reconstruct(const Vertices& g, const Vector& x) const
{
    validate_vertices(g);
    validate(x, blocks_->ncoarse(), "reconstruction x");
    Vector rhs;
    blocks_->apply_b(x, rhs);
    for (std::size_t i = 0; i < rhs.size(); ++i)
    {
        rhs[i] = -g.f[i] - rhs[i];
    }
    Vector result = solve_virtual(rhs);
    result.insert(result.end(), x.begin(), x.end());
    return result;
}

SternheimerWeakAugmented::Worker::SolveResult SternheimerWeakAugmented::Worker::solve(
    const Vertices& g, const SternheimerRPA::SolverOptions& options, const int restart_dimension) const
{
    validate_vertices(g);
    if (options.max_iter <= 0 || restart_dimension <= 0
        || !std::isfinite(options.residual_tol) || options.residual_tol <= 0.0
        || !std::isfinite(options.breakdown_tol) || options.breakdown_tol <= 0.0)
    {
        throw std::invalid_argument("weak augmented invalid GMRES options");
    }
    budget(workspace_bytes_required(restart_dimension), max_workspace_bytes_);
    SolveResult result;
    Vector x(blocks_->ncoarse(), 0.0);
    if (blocks_->ncoarse() != 0)
    {
        const Vector rhs = schur_rhs(g);
        Vector normalized_rhs;
        blocks_->apply_complement_inverse_sqrt(rhs, normalized_rhs);
        const double rhs_norm = norm(rhs);
        const double normalized_rhs_norm = norm(normalized_rhs);
        Vector source = g.f;
        source.insert(source.end(), g.w.begin(), g.w.end());
        auto normalized_options = options;
        if (normalized_rhs_norm > 0.0)
        {
            // S=I-L*L <= I implies ||r|| <= ||T r||. Target both original
            // Schur and full-source norms, with margin for the coordinate map.
            const double ratio = std::min(1.0, std::min(rhs_norm, norm(source)) / normalized_rhs_norm);
            normalized_options.residual_tol *= 0.5 * ratio;
            if (!(normalized_options.residual_tol > 0.0) || !std::isfinite(normalized_options.residual_tol))
                throw std::invalid_argument("weak augmented normalized residual target is not representable");
        }
        SternheimerRPA::LinearProblem problem;
        problem.apply = [this](const Vector& input, Vector& output) {
            Vector original, applied;
            blocks_->apply_complement_inverse_sqrt(input, original);
            apply_schur(original, applied);
            blocks_->apply_complement_inverse_sqrt(applied, output);
        };
        problem.dot = euclidean_dot;
        Vector normalized_x(blocks_->ncoarse(), 0.0);
        result.schur = SternheimerRPA::solve_gmres(problem, normalized_rhs, normalized_x, normalized_options, restart_dimension);
        blocks_->apply_complement_inverse_sqrt(normalized_x, x);
        // Existing GMRES may stop on its small least-squares residual. Always
        // reapply the actual operator before declaring even Schur convergence.
        Vector residual;
        apply_schur(x, residual);
        for (std::size_t i = 0; i < rhs.size(); ++i)
        {
            residual[i] -= rhs[i];
        }
        result.schur.absolute_residual = norm(residual);
        result.schur.relative_residual = result.schur.absolute_residual / (rhs_norm > 0.0 ? rhs_norm : 1.0);
        // The conservative transformed target may exhaust the iteration budget
        // after the original equation is already solved; the original gate wins.
        result.schur.converged = result.schur.relative_residual <= options.residual_tol;
    }
    else
    {
        result.schur.converged = true;
    }
    result.coefficients = reconstruct(g, x);
    Vector residual;
    apply_shifted(result.coefficients, residual);
    Vector source = g.f;
    source.insert(source.end(), g.w.begin(), g.w.end());
    for (std::size_t i = 0; i < source.size(); ++i)
    {
        residual[i] += source[i];
    }
    result.absolute_residual = norm(residual);
    const double source_norm = norm(source);
    result.relative_residual = result.absolute_residual / (source_norm > 0.0 ? source_norm : 1.0);
    result.converged = result.schur.converged && result.relative_residual <= options.residual_tol;
    return result;
}

} // namespace ModuleRI
