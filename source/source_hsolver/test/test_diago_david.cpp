// Public Davidson and Davidson-subspace stability regression tests.
// CTest runs one case per process and requires the final success marker: some
// historical failure paths returned exit(0) before producing valid eigenpairs.
#include "source_hsolver/diago_david.h"
#include "source_hsolver/diago_dav_subspace.h"
#include "source_hsolver/hs_operator.h"
#include "source_base/parallel_global.h"
#include "source_base/module_container/base/third_party/lapack.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
typedef std::complex<double> Complex;

const char* const case_names[] = {
    "exact", "partial_exact", "collinear", "near_collinear",
    "degenerate", "generalized_spd", "warm_start", "small_direction", "restart_spd"
};

struct Problem
{
    int dimension = 12;
    int bands = 2;
    int leading_dimension = 15;
    int subspace_multiplier = 2;
    int max_iterations = 200;
    int solves = 1;
    double solver_threshold = 1.0e-12;
    double eigenvalue_tolerance = 1.0e-10;
    double residual_tolerance = 1.0e-10;
    double gram_tolerance = 1.0e-10;
    std::vector<Complex> h;
    std::vector<Complex> s;
    std::vector<Complex> initial;
    std::vector<double> precondition;
};

void multiply(const std::vector<Complex>& matrix, const int dimension,
              const Complex* in, Complex* out, const int ld, const int columns)
{
    for (int column = 0; column < columns; ++column)
    {
        for (int row = 0; row < dimension; ++row)
        {
            Complex value(0.0, 0.0);
            for (int index = 0; index < dimension; ++index)
            {
                value += matrix[row + index * dimension] * in[index + column * ld];
            }
            out[row + column * ld] = value;
        }
    }
}

class DenseOperator : public hsolver::HSOperator<Complex, base_device::DEVICE_CPU>
{
  public:
    DenseOperator(const std::vector<Complex>& h, const std::vector<Complex>& s, const int dimension)
        : h_(h), s_(s), dimension_(dimension)
    {
    }

    void update_k(const int ik) override
    {
        if (ik != 0)
        {
            throw std::runtime_error("Only one algebraic k-point is available.");
        }
    }

    void hpsi(const Complex* in, Complex* out, const int ld, const int columns) const override
    {
        ++h_calls_;
        multiply(h_, dimension_, in, out, ld, columns);
    }

    int h_calls() const
    {
        return h_calls_;
    }

    void spsi(const Complex* in, Complex* out, const int ld, const int columns) const override
    {
        multiply(s_, dimension_, in, out, ld, columns);
    }

  private:
    const std::vector<Complex>& h_;
    const std::vector<Complex>& s_;
    const int dimension_;
    mutable int h_calls_ = 0;
};

void set_hermitian(Problem& problem, const int row, const int column, const Complex value)
{
    problem.h[row + column * problem.dimension] = value;
    problem.h[column + row * problem.dimension] = std::conj(value);
}

Problem make_problem(const std::string& name)
{
    Problem problem;
    if (name == "small_direction")
    {
        problem.bands = 1;
        problem.max_iterations = 50;
        problem.residual_tolerance = 1.0e-12;
        problem.eigenvalue_tolerance = 1.0e-12;
    }
    else if (name == "restart_spd")
    {
        problem.dimension = 24;
        problem.leading_dimension = 27;
        problem.bands = 3;
        problem.solver_threshold = 1.0e-14;
        problem.residual_tolerance = 2.0e-7;
    }
    const int matrix_size = problem.dimension * problem.dimension;
    const int block_size = problem.leading_dimension * problem.bands;
    const Complex zero(0.0, 0.0);
    problem.h.assign(matrix_size, zero);
    problem.s.assign(matrix_size, zero);
    problem.initial.assign(block_size, zero);
    problem.precondition.assign(problem.dimension, 1.0);
    for (int row = 0; row < problem.dimension; ++row)
    {
        problem.h[row + row * problem.dimension] = 10.0 + row;
        problem.s[row + row * problem.dimension] = 1.0;
    }
    problem.h[0] = 1.0;
    problem.h[1 + problem.dimension] = 2.0;
    problem.h[2 + 2 * problem.dimension] = 5.0;
    problem.initial[0] = 1.0;
    if (problem.bands > 1)
    {
        problem.initial[1 + problem.leading_dimension] = 1.0;
    }
    if (problem.bands > 2)
    {
        problem.initial[2 + 2 * problem.leading_dimension] = 1.0;
    }

    if (name == "exact")
    {
        // Both initial columns are exact lowest eigenvectors. Expanding a zero
        // residual must not cause normalization by zero or an abort.
    }
    else if (name == "partial_exact")
    {
        const Complex coupling(0.3, 0.0);
        set_hermitian(problem, 1, 2, coupling);
    }
    else if (name == "collinear" || name == "near_collinear" || name == "warm_start")
    {
        const Complex first_coupling(0.2, 0.0);
        const Complex second_coupling(0.3, 0.0);
        set_hermitian(problem, 0, 2, first_coupling);
        set_hermitian(problem, 1, 2, second_coupling);
        if (name == "near_collinear")
        {
            const Complex perturbation(1.0e-12, 0.0);
            set_hermitian(problem, 1, 3, perturbation);
        }
        if (name == "warm_start")
        {
            problem.solves = 2;
        }
    }
    else if (name == "degenerate")
    {
        problem.h[1 + problem.dimension] = 1.0;
        problem.h[3 + 3 * problem.dimension] = 6.0;
        const double first_norm = std::sqrt(1.04);
        const double second_norm = std::sqrt(1.09);
        problem.initial[0] = 1.0 / first_norm;
        problem.initial[2] = 0.2 / first_norm;
        problem.initial[1 + problem.leading_dimension] = 1.0 / second_norm;
        problem.initial[3 + problem.leading_dimension] = 0.3 / second_norm;
    }
    else if (name == "generalized_spd" || name == "restart_spd")
    {
        // H = sqrt(S) A sqrt(S), with non-real Hermitian couplings in A.
        // The two initial correction directions remain exactly dependent.
        const Complex first_coupling(0.0, 0.2);
        const Complex second_coupling(0.3, 0.1);
        set_hermitian(problem, 0, 2, first_coupling);
        set_hermitian(problem, 1, 2, second_coupling);
        if (name == "restart_spd")
        {
            // A dense Hermitian problem and a two-block budget require repeated
            // restarts. Independent LAPACK residual and S-Gram checks verify
            // that reuse of cached H/S products preserves the generalized solve.
            for (int column = 1; column < problem.dimension; ++column)
            {
                for (int row = 0; row < column; ++row)
                {
                    const double real_part = 0.15 / (column - row);
                    const double imaginary_part = 0.03 * (row + 1) / column;
                    const Complex coupling(real_part, imaginary_part);
                    set_hermitian(problem, row, column, coupling);
                }
            }
        }
        std::vector<double> scales(problem.dimension, 1.0);
        for (int row = 0; row < problem.dimension; ++row)
        {
            const double diagonal = 1.0 + 0.2 * row;
            scales[row] = std::sqrt(diagonal);
            problem.s[row + row * problem.dimension] = diagonal;
        }
        for (int column = 0; column < problem.dimension; ++column)
        {
            for (int row = 0; row < problem.dimension; ++row)
            {
                problem.h[row + column * problem.dimension] *= scales[row] * scales[column];
            }
        }
        for (int column = 0; column < problem.bands; ++column)
        {
            for (int row = 0; row < problem.dimension; ++row)
            {
                problem.initial[row + column * problem.leading_dimension] /= scales[row];
            }
        }
    }
    else if (name == "small_direction")
    {
        // True residual is about 2e-11, above ethr; preconditioning shrinks
        // the useful direction to about 2e-13. Absolute step-size rejection
        // must not accept the inaccurate input or abort before returning.
        for (int row = 0; row < problem.dimension; ++row)
        {
            problem.h[row + row * problem.dimension] = row;
        }
        problem.h[0] = -1.0;
        const double norm_squared = 1.0 + 1.0e-22;
        const double norm = std::sqrt(norm_squared);
        problem.initial[0] = 1.0 / norm;
        problem.initial[1] = 1.0e-11 / norm;
        problem.precondition.assign(problem.dimension, 100.0);
    }
    else
    {
        const std::string message = "Unknown case: " + name;
        throw std::runtime_error(message);
    }
    // Padding is deliberately nonzero: the operator contract uses only n rows.
    for (int column = 0; column < problem.bands; ++column)
    {
        for (int row = problem.dimension; row < problem.leading_dimension; ++row)
        {
            problem.initial[row + column * problem.leading_dimension] = Complex(17.0, -9.0);
        }
    }
    return problem;
}

struct Reference
{
    std::vector<double> eigenvalues;
    std::vector<Complex> vectors;
    int info = -999;
};

Reference reference_solve(const Problem& problem)
{
    Reference reference;
    const int n = problem.dimension;
    reference.eigenvalues.assign(n, 0.0);
    reference.vectors = problem.h;
    std::vector<Complex> overlap(problem.s);
    const int itype = 1;
    const char jobz = 'V';
    const char uplo = 'U';
    int lwork = -1;
    int lrwork = -1;
    int liwork = -1;
    Complex work_query(0.0, 0.0);
    double rwork_query = 0.0;
    int iwork_query = 0;
    zhegvd_(&itype, &jobz, &uplo, &n, reference.vectors.data(), &n, overlap.data(), &n,
            reference.eigenvalues.data(), &work_query, &lwork, &rwork_query, &lrwork,
            &iwork_query, &liwork, &reference.info);
    if (reference.info != 0)
    {
        return reference;
    }
    lwork = static_cast<int>(work_query.real());
    lrwork = static_cast<int>(rwork_query);
    liwork = iwork_query;
    if (lwork < 1 || lrwork < 1 || liwork < 1)
    {
        throw std::runtime_error("Invalid LAPACK workspace query result.");
    }
    std::vector<Complex> work(lwork);
    std::vector<double> rwork(lrwork);
    std::vector<int> iwork(liwork);
    reference.vectors = problem.h;
    overlap = problem.s;
    zhegvd_(&itype, &jobz, &uplo, &n, reference.vectors.data(), &n, overlap.data(), &n,
            reference.eigenvalues.data(), work.data(), &lwork, rwork.data(), &lrwork,
            iwork.data(), &liwork, &reference.info);
    return reference;
}

struct Metrics
{
    bool finite = true;
    double eigenvalue_error = 0.0;
    double residual = 0.0;
    double normalized_residual = 0.0;
    double gram_error = 0.0;
};

Metrics measure(const Problem& problem, const Complex* vectors, const int ld, const int columns,
                const std::vector<double>& eigenvalues, const std::vector<double>& reference)
{
    const int length = ld * columns;
    const Complex zero(0.0, 0.0);
    std::vector<Complex> hv(length, zero);
    std::vector<Complex> sv(length, zero);
    multiply(problem.h, problem.dimension, vectors, hv.data(), ld, columns);
    multiply(problem.s, problem.dimension, vectors, sv.data(), ld, columns);
    Metrics result;
    for (int column = 0; column < columns; ++column)
    {
        const double eigenvalue_difference = eigenvalues[column] - reference[column];
        const double delta = std::abs(eigenvalue_difference);
        result.finite = result.finite && std::isfinite(delta);
        result.eigenvalue_error = std::max(result.eigenvalue_error, delta);
        double residual2 = 0.0;
        double norm2 = 0.0;
        for (int row = 0; row < problem.dimension; ++row)
        {
            const int index = row + column * ld;
            const Complex difference = hv[index] - eigenvalues[column] * sv[index];
            residual2 += std::norm(difference);
            norm2 += std::norm(vectors[index]);
        }
        const double residual = std::sqrt(residual2);
        const double normalized_squared = residual2 / norm2;
        const double normalized = std::sqrt(normalized_squared);
        result.finite = result.finite && norm2 > 0.0 && std::isfinite(residual) && std::isfinite(normalized);
        result.residual = std::max(result.residual, residual);
        result.normalized_residual = std::max(result.normalized_residual, normalized);
    }
    for (int right = 0; right < columns; ++right)
    {
        for (int left = 0; left < columns; ++left)
        {
            Complex product(0.0, 0.0);
            for (int row = 0; row < problem.dimension; ++row)
            {
                const Complex left_value = vectors[row + left * ld];
                product += std::conj(left_value) * sv[row + right * ld];
            }
            if (left == right)
            {
                product -= 1.0;
            }
            const double error = std::abs(product);
            result.finite = result.finite && std::isfinite(error);
            result.gram_error = std::max(result.gram_error, error);
        }
    }
    return result;
}

void print_metrics(const std::string& marker, const std::string& solver, const std::string& name,
                   const int pass, const Metrics& metrics)
{
    std::cout << marker << " solver=" << solver << " case=" << name << " pass=" << pass
              << " finite=" << metrics.finite << " eigenvalue_error=" << metrics.eigenvalue_error
              << " residual=" << metrics.residual << " normalized_residual=" << metrics.normalized_residual
              << " SGram_error=" << metrics.gram_error << std::endl;
}

bool validate_return(const Problem& problem, const Reference& reference,
                     const std::vector<Complex>& guarded, const std::vector<double>& eigenvalues,
                     const std::string& solver, const std::string& name, const int pass,
                     const int iterations)
{
    const int guard = 3;
    const int block_size = problem.leading_dimension * problem.bands;
    const Complex sentinel(123.0, -456.0);
    bool canaries = true;
    for (int index = 0; index < guard; ++index)
    {
        canaries = canaries && guarded[index] == sentinel;
        canaries = canaries && guarded[guard + block_size + index] == sentinel;
    }
    const Complex* vectors = guarded.data() + guard;
    // Both public diag implementations explicitly clear the entire output
    // block before writing active rows. Their padding contract is zero, not
    // preservation of the nonzero input padding used by this test.
    bool padding = true;
    const Complex zero_padding(0.0, 0.0);
    for (int column = 0; column < problem.bands; ++column)
    {
        for (int row = problem.dimension; row < problem.leading_dimension; ++row)
        {
            const int index = row + column * problem.leading_dimension;
            padding = padding && vectors[index] == zero_padding;
        }
    }
    const Metrics metrics = measure(problem, vectors, problem.leading_dimension, problem.bands,
                                    eigenvalues, reference.eigenvalues);
    print_metrics("STABILITY_METRIC", solver, name, pass, metrics);
    bool ordered = true;
    for (int band = 0; band < problem.bands; ++band)
    {
        if (band > 0)
        {
            ordered = ordered && eigenvalues[band] >= eigenvalues[band - 1];
        }
        std::cout << "STABILITY_EIGENVALUE solver=" << solver << " case=" << name
                  << " pass=" << pass << " band=" << band + 1 << " value=" << eigenvalues[band]
                  << " reference=" << reference.eigenvalues[band] << std::endl;
    }
    const bool accepted = metrics.finite && canaries && padding && ordered && iterations >= 0
                          && metrics.eigenvalue_error <= problem.eigenvalue_tolerance
                          && metrics.residual <= problem.residual_tolerance
                          && metrics.gram_error <= problem.gram_tolerance;
    std::cout << "STABILITY_CHECK solver=" << solver << " case=" << name << " pass=" << pass
              << " accepted=" << accepted << " canaries=" << canaries
              << " padding=" << padding << " ordered=" << ordered
              << " iterations=" << iterations << std::endl;
    return accepted;
}

int run_case(const std::string& solver_name, const std::string& name,
             const hsolver::diag_comm_info& communicator)
{
    const Problem problem = make_problem(name);
    const Reference reference = reference_solve(problem);
    std::cout << "STABILITY_INPUT solver=" << solver_name << " case=" << name
              << " n=" << problem.dimension << " bands=" << problem.bands
              << " ld=" << problem.leading_dimension << " ndim=" << problem.subspace_multiplier
              << " ethr=" << problem.solver_threshold << " maxiter=" << problem.max_iterations
              << " solves=" << problem.solves << " eigenvalue_tolerance=" << problem.eigenvalue_tolerance
              << " residual_tolerance=" << problem.residual_tolerance
              << " gram_tolerance=" << problem.gram_tolerance << " reference_info=" << reference.info
              << " matrix_energy_unit=algebraic" << std::endl;
    if (reference.info != 0)
    {
        throw std::runtime_error("Independent direct LAPACK reference solve failed.");
    }
    const Metrics reference_metrics = measure(problem, reference.vectors.data(), problem.dimension,
                                              problem.dimension, reference.eigenvalues, reference.eigenvalues);
    print_metrics("STABILITY_REFERENCE", solver_name, name, 0, reference_metrics);
    if (!reference_metrics.finite || reference_metrics.residual > 1.0e-11 || reference_metrics.gram_error > 1.0e-11)
    {
        throw std::runtime_error("Independent direct LAPACK reference validation failed.");
    }
    const int guard = 3;
    const int block_size = problem.leading_dimension * problem.bands;
    const int guarded_size = block_size + 2 * guard;
    const Complex sentinel(123.0, -456.0);
    std::vector<Complex> guarded(guarded_size, sentinel);
    Complex* vectors = guarded.data() + guard;
    std::copy(problem.initial.begin(), problem.initial.end(), vectors);
    std::vector<double> eigenvalues(problem.bands, 0.0);
    const std::vector<double> thresholds(problem.bands, problem.solver_threshold);
    DenseOperator op(problem.h, problem.s, problem.dimension);
    op.update_k(0);
    bool accepted = true;
    if (solver_name == "dav")
    {
        hsolver::DiagoDavid<Complex, base_device::DEVICE_CPU> solver(problem.precondition.data(), problem.bands,
            problem.dimension, problem.subspace_multiplier, communicator);
        const int restart_budget = 5;
        const int allowed_unconverged = 0;
        for (int pass = 1; pass <= problem.solves; ++pass)
        {
            const int iterations = solver.diag(op, problem.leading_dimension, vectors, eigenvalues.data(),
                thresholds, problem.max_iterations, restart_budget, allowed_unconverged);
            std::cout << "STABILITY_RETURNED solver=" << solver_name << " case=" << name
                      << " pass=" << pass << std::endl;
            const bool pass_ok = validate_return(problem, reference, guarded, eigenvalues,
                                                 solver_name, name, pass, iterations);
            accepted = accepted && pass_ok;
        }
    }
    else if (solver_name == "dav_subspace")
    {
        const int lapack_path = 0;
        const int block_size_2d = 0;
        const bool scf_type = false;
        hsolver::Diago_DavSubspace<Complex, base_device::DEVICE_CPU> solver(problem.precondition, problem.bands,
            problem.dimension, problem.subspace_multiplier, problem.solver_threshold, problem.max_iterations,
            communicator, lapack_path, block_size_2d);
        for (int pass = 1; pass <= problem.solves; ++pass)
        {
            const int h_calls_before = op.h_calls();
            const int iterations = solver.diag(op, vectors, problem.leading_dimension, eigenvalues.data(),
                                               thresholds, scf_type);
            std::cout << "STABILITY_RETURNED solver=" << solver_name << " case=" << name
                      << " pass=" << pass << std::endl;
            const bool pass_ok = validate_return(problem, reference, guarded, eigenvalues,
                                                 solver_name, name, pass, iterations);
            accepted = accepted && pass_ok;
            if (name == "restart_spd")
            {
                const int h_calls = op.h_calls() - h_calls_before;
                // One initial H application plus at most one per expansion.
                // A restart of an accurate Ritz basis needs no extra H call.
                const bool reused_images = iterations >= 3 && h_calls <= iterations + 1;
                std::cout << "STABILITY_RESTART h_calls=" << h_calls << " iterations=" << iterations
                          << " reused_images=" << reused_images << std::endl;
                accepted = accepted && reused_images;
            }
        }
    }
    else
    {
        const std::string message = "Unknown solver: " + solver_name;
        throw std::runtime_error(message);
    }
    std::cout << "STABILITY_RESULT solver=" << solver_name << " case=" << name
              << " accepted=" << accepted << " completed_solves=" << problem.solves << std::endl;
    return accepted ? 0 : 1;
}
} // namespace

int main(int argc, char** argv)
{
    if (argc == 2 && std::string(argv[1]) == "--list")
    {
        for (const char* name : case_names)
        {
            std::cout << name << '\n';
        }
        return 0;
    }
    if (argc != 3)
    {
        std::cerr << "Usage: MODULE_HSOLVER_davidson_stability {dav|dav_subspace} CASE\n";
        return 2;
    }
    int nproc = 1;
    int rank = 0;
    int threads = 1;
    Parallel_Global::read_pal_param(argc, argv, nproc, threads, rank);
#ifdef __MPI
    // Existing reductions consume POOL_WORLD; the serial harness shares WORLD.
    // Other communicator handles are unused, and must be valid at finalize.
    POOL_WORLD = MPI_COMM_WORLD;
    KP_WORLD = MPI_COMM_NULL;
    BP_WORLD = MPI_COMM_NULL;
    INT_BGROUP = MPI_COMM_NULL;
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
    const hsolver::diag_comm_info communicator(POOL_WORLD, rank, nproc);
#else
    const hsolver::diag_comm_info communicator(rank, nproc);
#endif
    std::cout << std::scientific << std::setprecision(17);
    int status = 2;
    try
    {
        if (nproc != 1)
        {
            throw std::runtime_error("This algebraic harness requires exactly one MPI process.");
        }
        const std::string solver_name(argv[1]);
        const std::string case_name(argv[2]);
        status = run_case(solver_name, case_name, communicator);
    }
    catch (const std::exception& error)
    {
        std::cerr << "STABILITY_EXCEPTION " << error.what() << std::endl;
        status = 2;
    }
#ifdef __MPI
    Parallel_Global::finalize_mpi();
#endif
    return status;
}
