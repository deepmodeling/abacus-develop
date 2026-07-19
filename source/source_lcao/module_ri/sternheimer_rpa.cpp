#include "source_lcao/module_ri/sternheimer_rpa.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace ModuleRI
{

#ifdef __GREENX_MINIMAX
extern "C"
{
void gx_minimax_grid_frequency_wrp(int* num_points,
                                   double* e_min,
                                   double* e_max,
                                   double* omega_points,
                                   double* omega_weights,
                                   int* ierr);
}
#endif

namespace
{

constexpr std::int32_t kChi0V1Marker = -41073291;
constexpr std::int32_t kChi0V1ComplexFlag = 1;
constexpr std::int32_t kCoulombV1Marker = -20129433;

void assert_same_size(const SternheimerRPA::Vector& lhs, const SternheimerRPA::Vector& rhs, const char* context)
{
    if (lhs.size() != rhs.size())
    {
        throw std::invalid_argument(std::string(context) + " vector sizes differ.");
    }
}

double vector_norm(const SternheimerRPA::LinearProblem& problem, const SternheimerRPA::Vector& vec)
{
    const auto norm2 = problem.dot(vec, vec);
    return std::sqrt(std::max(0.0, norm2.real()));
}

void apply_preconditioner(const SternheimerRPA::LinearProblem& problem,
                          const SternheimerRPA::Vector& input,
                          SternheimerRPA::Vector& output)
{
    if (problem.precondition)
    {
        problem.precondition(input, output);
    }
    else
    {
        output = input;
    }
}

void axpy(const SternheimerRPA::Complex alpha, const SternheimerRPA::Vector& x, SternheimerRPA::Vector& y)
{
    assert_same_size(x, y, "SternheimerRPA::axpy");
    for (std::size_t i = 0; i != y.size(); ++i)
    {
        y[i] += alpha * x[i];
    }
}

void scale(const SternheimerRPA::Complex alpha, SternheimerRPA::Vector& x)
{
    for (SternheimerRPA::Complex& value: x)
    {
        value *= alpha;
    }
}

std::vector<SternheimerRPA::Complex> solve_small_dense_system(
    std::vector<std::vector<SternheimerRPA::Complex>> matrix,
    std::vector<SternheimerRPA::Complex> rhs)
{
    const int size = static_cast<int>(rhs.size());
    for (int pivot = 0; pivot != size; ++pivot)
    {
        int pivot_row = pivot;
        double pivot_abs = std::abs(matrix[pivot][pivot]);
        for (int row = pivot + 1; row != size; ++row)
        {
            const double candidate_abs = std::abs(matrix[row][pivot]);
            if (candidate_abs > pivot_abs)
            {
                pivot_abs = candidate_abs;
                pivot_row = row;
            }
        }
        if (pivot_abs < 1.0e-28)
        {
            throw std::runtime_error("SternheimerRPA::solve_gmres found a singular least-squares system.");
        }
        if (pivot_row != pivot)
        {
            std::swap(matrix[pivot], matrix[pivot_row]);
            std::swap(rhs[pivot], rhs[pivot_row]);
        }

        const SternheimerRPA::Complex diagonal = matrix[pivot][pivot];
        for (int col = pivot; col != size; ++col)
        {
            matrix[pivot][col] /= diagonal;
        }
        rhs[pivot] /= diagonal;

        for (int row = 0; row != size; ++row)
        {
            if (row == pivot)
            {
                continue;
            }
            const SternheimerRPA::Complex factor = matrix[row][pivot];
            if (std::abs(factor) == 0.0)
            {
                continue;
            }
            for (int col = pivot; col != size; ++col)
            {
                matrix[row][col] -= factor * matrix[pivot][col];
            }
            rhs[row] -= factor * rhs[pivot];
        }
    }
    return rhs;
}

struct LeastSquaresSolution
{
    std::vector<SternheimerRPA::Complex> coefficients;
    double residual_norm = 0.0;
};

LeastSquaresSolution solve_gmres_least_squares(
    const std::vector<std::vector<SternheimerRPA::Complex>>& hessenberg,
    const int rows,
    const int cols,
    const double beta)
{
    std::vector<std::vector<SternheimerRPA::Complex>> normal_matrix(
        cols, std::vector<SternheimerRPA::Complex>(cols, SternheimerRPA::Complex(0.0, 0.0)));
    std::vector<SternheimerRPA::Complex> normal_rhs(cols, SternheimerRPA::Complex(0.0, 0.0));

    for (int col = 0; col != cols; ++col)
    {
        normal_rhs[col] = std::conj(hessenberg[0][col]) * beta;
        for (int other = 0; other != cols; ++other)
        {
            for (int row = 0; row != rows; ++row)
            {
                normal_matrix[col][other] += std::conj(hessenberg[row][col]) * hessenberg[row][other];
            }
        }
    }

    LeastSquaresSolution result;
    result.coefficients = solve_small_dense_system(std::move(normal_matrix), std::move(normal_rhs));

    double residual2 = 0.0;
    for (int row = 0; row != rows; ++row)
    {
        SternheimerRPA::Complex value = row == 0 ? SternheimerRPA::Complex(beta, 0.0)
                                                 : SternheimerRPA::Complex(0.0, 0.0);
        for (int col = 0; col != cols; ++col)
        {
            value -= hessenberg[row][col] * result.coefficients[col];
        }
        residual2 += std::norm(value);
    }
    result.residual_norm = std::sqrt(residual2);
    return result;
}

void add_krylov_update(const std::vector<SternheimerRPA::Vector>& preconditioned_basis,
                       const std::vector<SternheimerRPA::Complex>& coefficients,
                       SternheimerRPA::Vector& solution)
{
    for (std::size_t basis_index = 0; basis_index != coefficients.size(); ++basis_index)
    {
        axpy(coefficients[basis_index], preconditioned_basis[basis_index], solution);
    }
}

void checked_write(std::ofstream& ofs, const void* data, const std::size_t bytes, const std::string& filename)
{
    const char* ptr = reinterpret_cast<const char*>(data);
    std::size_t bytes_left = bytes;
    const std::size_t max_chunk = static_cast<std::size_t>(std::numeric_limits<std::streamsize>::max());
    while (bytes_left > 0)
    {
        const std::size_t chunk = std::min(bytes_left, max_chunk);
        ofs.write(ptr, static_cast<std::streamsize>(chunk));
        if (!ofs.good())
        {
            throw std::runtime_error("Failed to write " + filename);
        }
        ptr += chunk;
        bytes_left -= chunk;
    }
}

template <typename Value>
void write_scalar(std::ofstream& ofs, const Value& value, const std::string& filename)
{
    checked_write(ofs, &value, sizeof(Value), filename);
}

void checked_read(std::ifstream& in, void* data, const std::size_t bytes, const std::string& filename)
{
    in.read(reinterpret_cast<char*>(data), static_cast<std::streamsize>(bytes));
    if (!in.good())
    {
        throw std::runtime_error("Failed to read " + filename);
    }
}

template <typename Value>
Value read_scalar(std::ifstream& in, const std::string& filename)
{
    Value value{};
    checked_read(in, &value, sizeof(Value), filename);
    return value;
}

std::int32_t checked_i32_from_size(const std::size_t value, const std::string& context)
{
    if (value > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
    {
        throw std::runtime_error(context + " exceeds int32_t range.");
    }
    return static_cast<std::int32_t>(value);
}

std::int64_t checked_i64_from_size(const std::size_t value, const std::string& context)
{
    if (value > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
    {
        throw std::runtime_error(context + " exceeds int64_t range.");
    }
    return static_cast<std::int64_t>(value);
}

std::size_t checked_mul_size(const std::size_t lhs, const std::size_t rhs, const std::string& context)
{
    if (lhs != 0 && rhs > std::numeric_limits<std::size_t>::max() / lhs)
    {
        throw std::runtime_error(context + " exceeds size_t range.");
    }
    return lhs * rhs;
}

int sum_positive_sizes(const std::vector<int>& values, const std::string& context)
{
    int sum = 0;
    for (const int value: values)
    {
        if (value <= 0)
        {
            throw std::invalid_argument(context + " contains a non-positive atom block size.");
        }
        if (value > std::numeric_limits<int>::max() - sum)
        {
            throw std::runtime_error(context + " exceeds int range.");
        }
        sum += value;
    }
    return sum;
}

std::size_t upper_triangular_pair_index(const std::size_t iatom, const std::size_t jatom, const std::size_t natom)
{
    if (iatom > jatom)
    {
        throw std::invalid_argument("Sternheimer chi0 v1 writer expects upper-triangular atom pairs.");
    }
    return iatom * natom - iatom * (iatom - 1) / 2 + (jatom - iatom);
}

std::pair<int, int> upper_triangular_pair_from_index(const int pair_index, const int natom)
{
    int index = 0;
    for (int iatom = 0; iatom != natom; ++iatom)
    {
        for (int jatom = iatom; jatom != natom; ++jatom)
        {
            if (index == pair_index)
            {
                return {iatom, jatom};
            }
            ++index;
        }
    }
    throw std::runtime_error("Coulomb v1 atom-pair index is out of range.");
}

std::vector<std::vector<int>> map_atom_local_to_global(
    const std::vector<SternheimerRPA::AuxiliaryChannel>& channels,
    const std::vector<int>& atom_naux)
{
    std::vector<std::vector<int>> global_index(atom_naux.size());
    for (std::size_t iatom = 0; iatom != atom_naux.size(); ++iatom)
    {
        global_index[iatom].assign(static_cast<std::size_t>(atom_naux[iatom]), -1);
    }

    for (std::size_t iglobal = 0; iglobal != channels.size(); ++iglobal)
    {
        const SternheimerRPA::AuxiliaryChannel& channel = channels[iglobal];
        if (channel.channel_index != static_cast<int>(iglobal))
        {
            throw std::invalid_argument("Sternheimer chi0 v1 channel indices must match matrix ordering.");
        }
        if (channel.atom_index < 0 || channel.atom_index >= static_cast<int>(atom_naux.size()))
        {
            throw std::invalid_argument("Sternheimer chi0 v1 channel atom index is out of range.");
        }
        const int naux = atom_naux[static_cast<std::size_t>(channel.atom_index)];
        if (channel.atom_local_index < 0 || channel.atom_local_index >= naux)
        {
            throw std::invalid_argument("Sternheimer chi0 v1 channel local auxiliary index is out of range.");
        }
        int& mapped_index
            = global_index[static_cast<std::size_t>(channel.atom_index)][static_cast<std::size_t>(
                channel.atom_local_index)];
        if (mapped_index != -1)
        {
            throw std::invalid_argument("Sternheimer chi0 v1 channel mapping contains a duplicate local index.");
        }
        mapped_index = static_cast<int>(iglobal);
    }

    for (const std::vector<int>& atom_map: global_index)
    {
        for (const int mapped_index: atom_map)
        {
            if (mapped_index < 0)
            {
                throw std::invalid_argument("Sternheimer chi0 v1 channel mapping is incomplete.");
            }
        }
    }
    return global_index;
}

struct Chi0Block
{
    int pair_index = 0;
    std::int64_t offset = 0;
    std::vector<SternheimerRPA::Complex> payload;
};

struct CoulombFileBlock
{
    int pair_index = 0;
    std::int64_t offset = 0;
};

} // namespace

SternheimerRPA::CoulombV1Matrix SternheimerRPA::read_coulomb_v1_files(
    const std::vector<std::string>& filenames)
{
    if (filenames.empty())
    {
        throw std::invalid_argument("Coulomb v1 reader requires at least one rank file.");
    }

    CoulombV1Matrix result;
    int value_flag_reference = -1;
    std::vector<bool> pair_seen;
    for (const std::string& filename: filenames)
    {
        std::ifstream in(filename.c_str(), std::ios::binary);
        if (!in.good())
        {
            throw std::runtime_error("Failed to open " + filename);
        }
        const std::int32_t marker = read_scalar<std::int32_t>(in, filename);
        const std::int32_t iq = read_scalar<std::int32_t>(in, filename);
        const std::int32_t naux = read_scalar<std::int32_t>(in, filename);
        const std::int32_t value_flag = read_scalar<std::int32_t>(in, filename);
        const std::int32_t natom = read_scalar<std::int32_t>(in, filename);
        const std::int32_t nblocks = read_scalar<std::int32_t>(in, filename);
        if (marker != kCoulombV1Marker || iq <= 0 || naux <= 0 || natom <= 0 || nblocks < 0)
        {
            throw std::runtime_error("Invalid Coulomb v1 header in " + filename);
        }
        std::vector<int> atom_naux(static_cast<std::size_t>(natom));
        for (int& count: atom_naux)
        {
            count = read_scalar<std::int32_t>(in, filename);
        }
        if (sum_positive_sizes(atom_naux, "Coulomb v1 atom sizes") != naux)
        {
            throw std::runtime_error("Coulomb v1 atom sizes do not sum to naux.");
        }
        if (result.iq == 0)
        {
            result.iq = iq;
            result.atom_naux = atom_naux;
            result.values.assign(static_cast<std::size_t>(naux) * static_cast<std::size_t>(naux), Complex(0.0, 0.0));
            value_flag_reference = value_flag;
            pair_seen.assign(static_cast<std::size_t>(natom * (natom + 1) / 2), false);
        }
        else if (result.iq != iq || result.atom_naux != atom_naux || value_flag_reference != value_flag)
        {
            throw std::runtime_error("Coulomb v1 rank-file metadata are inconsistent.");
        }

        std::vector<CoulombFileBlock> blocks(static_cast<std::size_t>(nblocks));
        for (CoulombFileBlock& block: blocks)
        {
            block.pair_index = read_scalar<std::int32_t>(in, filename);
            block.offset = read_scalar<std::int64_t>(in, filename);
        }
        std::vector<int> starts(static_cast<std::size_t>(natom), 0);
        for (int iatom = 1; iatom != natom; ++iatom)
        {
            starts[static_cast<std::size_t>(iatom)]
                = starts[static_cast<std::size_t>(iatom - 1)] + atom_naux[static_cast<std::size_t>(iatom - 1)];
        }
        for (const CoulombFileBlock& block: blocks)
        {
            if (block.pair_index < 0 || block.pair_index >= static_cast<int>(pair_seen.size())
                || pair_seen[static_cast<std::size_t>(block.pair_index)])
            {
                throw std::runtime_error("Coulomb v1 atom-pair blocks are invalid or duplicated.");
            }
            pair_seen[static_cast<std::size_t>(block.pair_index)] = true;
            const auto pair = upper_triangular_pair_from_index(block.pair_index, natom);
            const int nrow = atom_naux[static_cast<std::size_t>(pair.first)];
            const int ncol = atom_naux[static_cast<std::size_t>(pair.second)];
            in.seekg(block.offset, std::ios::beg);
            if (!in.good())
            {
                throw std::runtime_error("Invalid Coulomb v1 block offset in " + filename);
            }
            for (int irow = 0; irow != nrow; ++irow)
            {
                for (int icol = 0; icol != ncol; ++icol)
                {
                    Complex value;
                    if (value_flag == 1)
                    {
                        value = read_scalar<Complex>(in, filename);
                    }
                    else if (value_flag == 0)
                    {
                        value = Complex(read_scalar<double>(in, filename), 0.0);
                    }
                    else
                    {
                        throw std::runtime_error("Unsupported Coulomb v1 value flag.");
                    }
                    const int row = starts[static_cast<std::size_t>(pair.first)] + irow;
                    const int col = starts[static_cast<std::size_t>(pair.second)] + icol;
                    result.values[static_cast<std::size_t>(row) * static_cast<std::size_t>(naux)
                                  + static_cast<std::size_t>(col)] = value;
                    if (pair.first != pair.second)
                    {
                        result.values[static_cast<std::size_t>(col) * static_cast<std::size_t>(naux)
                                      + static_cast<std::size_t>(row)] = std::conj(value);
                    }
                }
            }
        }
    }
    if (std::find(pair_seen.begin(), pair_seen.end(), false) != pair_seen.end())
    {
        throw std::runtime_error("Coulomb v1 rank files do not contain every atom-pair block.");
    }
    return result;
}

SternheimerRPA::SolverResult SternheimerRPA::solve_bicgstab(const LinearProblem& problem,
                                                            const Vector& rhs,
                                                            Vector& solution)
{
    return solve_bicgstab(problem, rhs, solution, SolverOptions{});
}

SternheimerRPA::SolverResult SternheimerRPA::solve_bicgstab(const LinearProblem& problem,
                                                            const Vector& rhs,
                                                            Vector& solution,
                                                            const SolverOptions& options)
{
    if (!problem.apply)
    {
        throw std::invalid_argument("SternheimerRPA::solve_bicgstab requires an apply callback.");
    }
    if (!problem.dot)
    {
        throw std::invalid_argument("SternheimerRPA::solve_bicgstab requires a dot callback.");
    }
    if (solution.empty())
    {
        solution.assign(rhs.size(), Complex(0.0, 0.0));
    }
    assert_same_size(rhs, solution, "SternheimerRPA::solve_bicgstab");

    Vector residual(rhs.size());
    Vector ax(rhs.size());
    problem.apply(solution, ax);
    for (std::size_t i = 0; i != rhs.size(); ++i)
    {
        residual[i] = rhs[i] - ax[i];
    }

    const double rhs_norm = vector_norm(problem, rhs);
    const double norm_floor = rhs_norm > 0.0 ? rhs_norm : 1.0;

    SolverResult result;
    result.absolute_residual = vector_norm(problem, residual);
    result.relative_residual = result.absolute_residual / norm_floor;
    if (result.relative_residual <= options.residual_tol)
    {
        result.converged = true;
        return result;
    }

    Vector residual_shadow = residual;
    Vector p(rhs.size(), Complex(0.0, 0.0));
    Vector v(rhs.size(), Complex(0.0, 0.0));
    Vector s(rhs.size());
    Vector t(rhs.size());
    Vector p_hat(rhs.size());
    Vector s_hat(rhs.size());

    Complex rho_old(1.0, 0.0);
    Complex alpha(1.0, 0.0);
    Complex omega_old(1.0, 0.0);

    for (int iter = 1; iter <= options.max_iter; ++iter)
    {
        const Complex rho_new = problem.dot(residual_shadow, residual);
        if (std::abs(rho_new) < options.breakdown_tol)
        {
            result.iterations = iter - 1;
            return result;
        }

        if (iter == 1)
        {
            p = residual;
        }
        else
        {
            const Complex beta = (rho_new / rho_old) * (alpha / omega_old);
            for (std::size_t i = 0; i != p.size(); ++i)
            {
                p[i] = residual[i] + beta * (p[i] - omega_old * v[i]);
            }
        }

        apply_preconditioner(problem, p, p_hat);
        problem.apply(p_hat, v);
        const Complex rv = problem.dot(residual_shadow, v);
        if (std::abs(rv) < options.breakdown_tol)
        {
            result.iterations = iter - 1;
            return result;
        }
        alpha = rho_new / rv;

        for (std::size_t i = 0; i != s.size(); ++i)
        {
            s[i] = residual[i] - alpha * v[i];
        }
        const double s_norm = vector_norm(problem, s);
        if (s_norm / norm_floor <= options.residual_tol)
        {
            axpy(alpha, p_hat, solution);
            result.converged = true;
            result.iterations = iter;
            result.absolute_residual = s_norm;
            result.relative_residual = s_norm / norm_floor;
            return result;
        }

        apply_preconditioner(problem, s, s_hat);
        problem.apply(s_hat, t);
        const Complex tt = problem.dot(t, t);
        if (std::abs(tt) < options.breakdown_tol)
        {
            result.iterations = iter - 1;
            return result;
        }
        const Complex omega_new = problem.dot(t, s) / tt;
        if (std::abs(omega_new) < options.breakdown_tol)
        {
            result.iterations = iter - 1;
            return result;
        }

        axpy(alpha, p_hat, solution);
        axpy(omega_new, s_hat, solution);

        for (std::size_t i = 0; i != residual.size(); ++i)
        {
            residual[i] = s[i] - omega_new * t[i];
        }

        result.iterations = iter;
        result.absolute_residual = vector_norm(problem, residual);
        result.relative_residual = result.absolute_residual / norm_floor;
        if (result.relative_residual <= options.residual_tol)
        {
            result.converged = true;
            return result;
        }

        rho_old = rho_new;
        omega_old = omega_new;
    }

    return result;
}

SternheimerRPA::SolverResult SternheimerRPA::solve_gmres(const LinearProblem& problem,
                                                         const Vector& rhs,
                                                         Vector& solution,
                                                         const SolverOptions& options,
                                                         const int restart_dimension)
{
    if (!problem.apply)
    {
        throw std::invalid_argument("SternheimerRPA::solve_gmres requires an apply callback.");
    }
    if (!problem.dot)
    {
        throw std::invalid_argument("SternheimerRPA::solve_gmres requires a dot callback.");
    }
    if (restart_dimension <= 0)
    {
        throw std::invalid_argument("SternheimerRPA::solve_gmres requires a positive restart dimension.");
    }
    if (solution.empty())
    {
        solution.assign(rhs.size(), Complex(0.0, 0.0));
    }
    assert_same_size(rhs, solution, "SternheimerRPA::solve_gmres");

    SolverResult result;
    Vector residual(rhs.size());
    Vector ax(rhs.size());
    problem.apply(solution, ax);
    for (std::size_t i = 0; i != rhs.size(); ++i)
    {
        residual[i] = rhs[i] - ax[i];
    }

    const double rhs_norm = vector_norm(problem, rhs);
    const double norm_floor = rhs_norm > 0.0 ? rhs_norm : 1.0;
    double beta = vector_norm(problem, residual);
    result.absolute_residual = beta;
    result.relative_residual = beta / norm_floor;
    if (result.relative_residual <= options.residual_tol)
    {
        result.converged = true;
        return result;
    }

    while (result.iterations < options.max_iter)
    {
        const int max_inner = std::min(restart_dimension, options.max_iter - result.iterations);
        std::vector<Vector> basis(max_inner + 1, Vector(rhs.size(), Complex(0.0, 0.0)));
        std::vector<Vector> preconditioned_basis(max_inner, Vector(rhs.size(), Complex(0.0, 0.0)));
        std::vector<std::vector<Complex>> hessenberg(max_inner + 1,
                                                     std::vector<Complex>(max_inner, Complex(0.0, 0.0)));

        basis[0] = residual;
        scale(Complex(1.0 / beta, 0.0), basis[0]);

        LeastSquaresSolution least_squares;
        int used_inner = 0;
        for (int inner = 0; inner != max_inner; ++inner)
        {
            apply_preconditioner(problem, basis[inner], preconditioned_basis[inner]);

            Vector work;
            problem.apply(preconditioned_basis[inner], work);
            for (int basis_index = 0; basis_index <= inner; ++basis_index)
            {
                hessenberg[basis_index][inner] = problem.dot(basis[basis_index], work);
                axpy(-hessenberg[basis_index][inner], basis[basis_index], work);
            }

            const double next_norm = vector_norm(problem, work);
            hessenberg[inner + 1][inner] = Complex(next_norm, 0.0);
            if (next_norm > options.breakdown_tol)
            {
                basis[inner + 1] = work;
                scale(Complex(1.0 / next_norm, 0.0), basis[inner + 1]);
            }

            used_inner = inner + 1;
            least_squares = solve_gmres_least_squares(hessenberg, used_inner + 1, used_inner, beta);
            result.absolute_residual = least_squares.residual_norm;
            result.relative_residual = result.absolute_residual / norm_floor;
            result.iterations += 1;

            if (result.relative_residual <= options.residual_tol || next_norm <= options.breakdown_tol)
            {
                add_krylov_update(preconditioned_basis, least_squares.coefficients, solution);
                result.converged = result.relative_residual <= options.residual_tol;
                return result;
            }
        }

        add_krylov_update(preconditioned_basis, least_squares.coefficients, solution);
        problem.apply(solution, ax);
        for (std::size_t i = 0; i != rhs.size(); ++i)
        {
            residual[i] = rhs[i] - ax[i];
        }
        beta = vector_norm(problem, residual);
        result.absolute_residual = beta;
        result.relative_residual = result.absolute_residual / norm_floor;
        if (result.relative_residual <= options.residual_tol)
        {
            result.converged = true;
            return result;
        }
    }

    return result;
}

void SternheimerRPA::build_rhs_from_hartree_perturbation(const std::vector<double>& hartree_potential_r,
                                                         const Vector& psi_r,
                                                         Vector& rhs_r)
{
    if (hartree_potential_r.size() != psi_r.size())
    {
        throw std::invalid_argument("SternheimerRPA::build_rhs_from_hartree_perturbation size mismatch.");
    }
    rhs_r.resize(psi_r.size());
    for (std::size_t ir = 0; ir != psi_r.size(); ++ir)
    {
        rhs_r[ir] = -hartree_potential_r[ir] * psi_r[ir];
    }
}

void SternheimerRPA::build_rhs_from_hartree_perturbation(const Vector& hartree_potential_r,
                                                         const Vector& psi_r,
                                                         Vector& rhs_r)
{
    if (hartree_potential_r.size() != psi_r.size())
    {
        throw std::invalid_argument("SternheimerRPA::build_rhs_from_hartree_perturbation size mismatch.");
    }
    rhs_r.resize(psi_r.size());
    for (std::size_t ir = 0; ir != psi_r.size(); ++ir)
    {
        rhs_r[ir] = -hartree_potential_r[ir] * psi_r[ir];
    }
}

SternheimerRPA::Complex SternheimerRPA::accumulate_polarizability_grid_element(
    const std::vector<double>& hartree_potential_r,
    const Vector& psi_r,
    const Vector& delta_psi_r,
    const double grid_weight)
{
    if (hartree_potential_r.size() != psi_r.size() || psi_r.size() != delta_psi_r.size())
    {
        throw std::invalid_argument("SternheimerRPA::accumulate_polarizability_grid_element size mismatch.");
    }
    Complex value(0.0, 0.0);
    for (std::size_t ir = 0; ir != psi_r.size(); ++ir)
    {
        value += std::conj(psi_r[ir]) * hartree_potential_r[ir] * delta_psi_r[ir];
    }
    return grid_weight * value;
}

SternheimerRPA::Complex SternheimerRPA::accumulate_polarizability_grid_element(
    const Vector& hartree_potential_r,
    const Vector& psi_r,
    const Vector& delta_psi_r,
    const double grid_weight)
{
    if (hartree_potential_r.size() != psi_r.size() || psi_r.size() != delta_psi_r.size())
    {
        throw std::invalid_argument("SternheimerRPA::accumulate_polarizability_grid_element size mismatch.");
    }
    Complex value(0.0, 0.0);
    for (std::size_t ir = 0; ir != psi_r.size(); ++ir)
    {
        value += std::conj(psi_r[ir]) * std::conj(hartree_potential_r[ir]) * delta_psi_r[ir];
    }
    return grid_weight * value;
}

void SternheimerRPA::accumulate_chi0_branch_column(const std::vector<std::vector<double>>& hartree_potentials_r,
                                                   const Vector& psi_r,
                                                   const Vector& delta_psi_r,
                                                   const double grid_weight,
                                                   const double occupation,
                                                   const int column_index,
                                                   std::vector<Complex>& branch_matrix)
{
    const int num_channels = static_cast<int>(hartree_potentials_r.size());
    if (num_channels <= 0)
    {
        throw std::invalid_argument("SternheimerRPA::accumulate_chi0_branch_column requires at least one channel.");
    }
    if (column_index < 0 || column_index >= num_channels)
    {
        throw std::invalid_argument("SternheimerRPA::accumulate_chi0_branch_column column index is out of range.");
    }
    if (branch_matrix.size()
        != checked_mul_size(static_cast<std::size_t>(num_channels),
                            static_cast<std::size_t>(num_channels),
                            "Sternheimer chi0 branch matrix size"))
    {
        throw std::invalid_argument("SternheimerRPA::accumulate_chi0_branch_column matrix size mismatch.");
    }

    for (int row = 0; row != num_channels; ++row)
    {
        branch_matrix[static_cast<std::size_t>(row) * static_cast<std::size_t>(num_channels)
                      + static_cast<std::size_t>(column_index)]
            += occupation * accumulate_polarizability_grid_element(hartree_potentials_r[static_cast<std::size_t>(row)],
                                                                   psi_r,
                                                                   delta_psi_r,
                                                                   grid_weight);
    }
}

void SternheimerRPA::accumulate_chi0_branch_column(const std::vector<Vector>& hartree_potentials_r,
                                                   const Vector& psi_r,
                                                   const Vector& delta_psi_r,
                                                   const double grid_weight,
                                                   const double occupation,
                                                   const int column_index,
                                                   std::vector<Complex>& branch_matrix)
{
    const int num_channels = static_cast<int>(hartree_potentials_r.size());
    if (num_channels <= 0)
    {
        throw std::invalid_argument("SternheimerRPA::accumulate_chi0_branch_column requires at least one channel.");
    }
    if (column_index < 0 || column_index >= num_channels)
    {
        throw std::invalid_argument("SternheimerRPA::accumulate_chi0_branch_column column index is out of range.");
    }
    if (branch_matrix.size()
        != checked_mul_size(static_cast<std::size_t>(num_channels),
                            static_cast<std::size_t>(num_channels),
                            "Sternheimer chi0 branch matrix size"))
    {
        throw std::invalid_argument("SternheimerRPA::accumulate_chi0_branch_column matrix size mismatch.");
    }

    for (int row = 0; row != num_channels; ++row)
    {
        branch_matrix[static_cast<std::size_t>(row) * static_cast<std::size_t>(num_channels)
                      + static_cast<std::size_t>(column_index)]
            += occupation * accumulate_polarizability_grid_element(hartree_potentials_r[static_cast<std::size_t>(row)],
                                                                   psi_r,
                                                                   delta_psi_r,
                                                                   grid_weight);
    }
}

std::vector<SternheimerRPA::Complex> SternheimerRPA::symmetrize_chi0_imaginary_frequency(
    const std::vector<Complex>& branch_matrix,
    const int num_channels)
{
    if (num_channels <= 0)
    {
        throw std::invalid_argument("SternheimerRPA::symmetrize_chi0_imaginary_frequency requires channels.");
    }
    if (branch_matrix.size()
        != checked_mul_size(static_cast<std::size_t>(num_channels),
                            static_cast<std::size_t>(num_channels),
                            "Sternheimer chi0 branch matrix size"))
    {
        throw std::invalid_argument("SternheimerRPA::symmetrize_chi0_imaginary_frequency matrix size mismatch.");
    }

    std::vector<Complex> chi0(branch_matrix.size(), Complex(0.0, 0.0));
    for (int row = 0; row != num_channels; ++row)
    {
        for (int col = 0; col != num_channels; ++col)
        {
            const std::size_t index
                = static_cast<std::size_t>(row) * static_cast<std::size_t>(num_channels)
                  + static_cast<std::size_t>(col);
            const std::size_t transpose
                = static_cast<std::size_t>(col) * static_cast<std::size_t>(num_channels)
                  + static_cast<std::size_t>(row);
            chi0[index] = branch_matrix[index] + std::conj(branch_matrix[transpose]);
        }
    }
    return chi0;
}

std::int32_t SternheimerRPA::chi0_v1_marker()
{
    return kChi0V1Marker;
}

int SternheimerRPA::frequency_owner_rank(const int ifrequency_zero_based,
                                         const int mpi_ranks,
                                         const int rank_shift)
{
    if (ifrequency_zero_based < 0)
    {
        throw std::invalid_argument("SternheimerRPA::frequency_owner_rank requires a non-negative frequency index.");
    }
    if (mpi_ranks <= 0)
    {
        throw std::invalid_argument("SternheimerRPA::frequency_owner_rank requires a positive MPI rank count.");
    }
    if (mpi_ranks == 1)
    {
        return 0;
    }

    int normalized_shift = rank_shift % mpi_ranks;
    if (normalized_shift < 0)
    {
        normalized_shift += mpi_ranks;
    }
    return (ifrequency_zero_based + normalized_shift) % mpi_ranks;
}

SternheimerRPA::TransitionEnergyWindow SternheimerRPA::transition_energy_window_from_eigenvalues_ry(
    const std::vector<double>& eigenvalues_ry,
    const std::vector<double>& occupations,
    const double occupation_tolerance)
{
    if (eigenvalues_ry.size() != occupations.size())
    {
        throw std::invalid_argument("Sternheimer minimax window requires matching eigenvalue and occupation sizes.");
    }
    if (eigenvalues_ry.empty())
    {
        throw std::invalid_argument("Sternheimer minimax window requires at least one eigenvalue.");
    }
    if (occupation_tolerance < 0.0)
    {
        throw std::invalid_argument("Sternheimer minimax occupation tolerance must be non-negative.");
    }

    double emin_ry = std::numeric_limits<double>::max();
    double emax_ry = 0.0;
    bool found_transition = false;
    for (std::size_t iv = 0; iv != eigenvalues_ry.size(); ++iv)
    {
        if (occupations[iv] <= occupation_tolerance)
        {
            continue;
        }
        for (std::size_t ic = 0; ic != eigenvalues_ry.size(); ++ic)
        {
            if (occupations[ic] > occupation_tolerance)
            {
                continue;
            }
            const double transition_ry = eigenvalues_ry[ic] - eigenvalues_ry[iv];
            if (transition_ry <= 0.0)
            {
                continue;
            }
            emin_ry = std::min(emin_ry, transition_ry);
            emax_ry = std::max(emax_ry, transition_ry);
            found_transition = true;
        }
    }

    if (!found_transition)
    {
        throw std::runtime_error("Sternheimer minimax window found no positive occupied-to-empty transition.");
    }

    TransitionEnergyWindow window;
    window.emin_ha = 0.5 * emin_ry;
    window.emax_ha = 0.5 * emax_ry;
    return window;
}

SternheimerRPA::FrequencyGrid SternheimerRPA::generate_greenx_minimax_frequency_grid(const int nfreq,
                                                                                     const double emin_ha,
                                                                                     const double emax_ha)
{
    if (nfreq <= 0)
    {
        throw std::invalid_argument("Sternheimer minimax grid requires nfreq > 0.");
    }
    if (emin_ha <= 0.0 || emax_ha <= 0.0 || emax_ha < emin_ha)
    {
        throw std::invalid_argument("Sternheimer minimax grid requires 0 < emin <= emax in Hartree.");
    }

    FrequencyGrid grid;
    grid.omega_ha.assign(static_cast<std::size_t>(nfreq), 0.0);
    grid.weights_ha.assign(static_cast<std::size_t>(nfreq), 0.0);

#ifdef __GREENX_MINIMAX
    int nfreq_c = nfreq;
    double emin_c = emin_ha;
    double emax_c = emax_ha;
    int ierr = -1;
    gx_minimax_grid_frequency_wrp(&nfreq_c,
                                  &emin_c,
                                  &emax_c,
                                  grid.omega_ha.data(),
                                  grid.weights_ha.data(),
                                  &ierr);
    if (ierr != 0)
    {
        throw std::runtime_error("GreenX minimax frequency-grid generation failed with ierr="
                                 + std::to_string(ierr) + ".");
    }
#else
    throw std::runtime_error("ABACUS was built without ENABLE_GREENX_MINIMAX; Sternheimer minimax is unavailable.");
#endif

    return grid;
}

SternheimerRPA::FrequencyGrid SternheimerRPA::read_frequency_grid_file(const std::string& filename,
                                                                       const int expected_size)
{
    if (expected_size <= 0)
    {
        throw std::invalid_argument("Sternheimer frequency grid file requires expected_size > 0.");
    }

    std::ifstream in(filename.c_str());
    if (!in.good())
    {
        throw std::runtime_error("Failed to open Sternheimer frequency grid file " + filename + ".");
    }

    FrequencyGrid grid;
    std::string line;
    int line_number = 0;
    while (std::getline(in, line))
    {
        ++line_number;
        const std::size_t comment = line.find('#');
        if (comment != std::string::npos)
        {
            line.resize(comment);
        }

        std::istringstream iss(line);
        std::vector<double> columns;
        double value = 0.0;
        while (iss >> value)
        {
            columns.push_back(value);
        }
        if (!iss.eof())
        {
            throw std::invalid_argument("Invalid token in Sternheimer frequency grid file " + filename
                                        + " at line " + std::to_string(line_number) + ".");
        }
        if (columns.empty())
        {
            continue;
        }
        if (columns.size() != 2 && columns.size() != 3)
        {
            throw std::invalid_argument("Sternheimer frequency grid file " + filename
                                        + " expects rows with 2 columns (omega_Ha weight_Ha) or 3 columns "
                                          "(index omega_Ha weight_Ha).");
        }

        double omega_ha = columns[0];
        double weight_ha = columns[1];
        if (columns.size() == 3)
        {
            const int expected_index = static_cast<int>(grid.omega_ha.size()) + 1;
            if (std::abs(columns[0] - static_cast<double>(expected_index)) > 1.0e-10)
            {
                throw std::invalid_argument("Sternheimer frequency grid file " + filename
                                            + " has an unexpected one-based index at line "
                                            + std::to_string(line_number) + ".");
            }
            omega_ha = columns[1];
            weight_ha = columns[2];
        }
        if (omega_ha <= 0.0 || weight_ha <= 0.0)
        {
            throw std::invalid_argument("Sternheimer frequency grid file " + filename
                                        + " requires positive omega_Ha and weight_Ha.");
        }

        grid.omega_ha.push_back(omega_ha);
        grid.weights_ha.push_back(weight_ha);
    }

    if (static_cast<int>(grid.omega_ha.size()) != expected_size)
    {
        throw std::runtime_error("Sternheimer frequency grid file " + filename + " contains "
                                 + std::to_string(grid.omega_ha.size()) + " rows but sternheimer_nfreq is "
                                 + std::to_string(expected_size) + ".");
    }
    return grid;
}

void SternheimerRPA::write_chi0_v1_file(const std::string& filename,
                                        const Chi0V1Metadata& metadata,
                                        const std::vector<AuxiliaryChannel>& channels,
                                        const std::vector<Complex>& chi0_matrix)
{
    if (metadata.iq <= 0 || metadata.ifrequency <= 0)
    {
        throw std::invalid_argument("Sternheimer chi0 v1 metadata requires positive iq and frequency indices.");
    }
    if (metadata.atom_naux.empty())
    {
        throw std::invalid_argument("Sternheimer chi0 v1 metadata requires atom auxiliary sizes.");
    }
    const int naux = sum_positive_sizes(metadata.atom_naux, "Sternheimer chi0 v1 atom_naux");
    if (channels.size() != static_cast<std::size_t>(naux))
    {
        throw std::invalid_argument("Sternheimer chi0 v1 channel count does not match atom_naux.");
    }
    const std::size_t matrix_size = checked_mul_size(static_cast<std::size_t>(naux),
                                                     static_cast<std::size_t>(naux),
                                                     "Sternheimer chi0 v1 matrix size");
    if (chi0_matrix.size() != matrix_size)
    {
        throw std::invalid_argument("Sternheimer chi0 v1 matrix size mismatch.");
    }

    const std::vector<std::vector<int>> atom_local_to_global
        = map_atom_local_to_global(channels, metadata.atom_naux);
    const std::size_t natom = metadata.atom_naux.size();
    std::vector<Chi0Block> blocks;
    blocks.reserve(natom * (natom + 1) / 2);
    for (std::size_t iatom = 0; iatom != natom; ++iatom)
    {
        for (std::size_t jatom = iatom; jatom != natom; ++jatom)
        {
            Chi0Block block;
            block.pair_index = checked_i32_from_size(upper_triangular_pair_index(iatom, jatom, natom),
                                                     "Sternheimer chi0 v1 atom-pair index");
            const int inaux = metadata.atom_naux[iatom];
            const int jnaux = metadata.atom_naux[jatom];
            block.payload.reserve(static_cast<std::size_t>(inaux) * static_cast<std::size_t>(jnaux));
            for (int imu = 0; imu != inaux; ++imu)
            {
                const int iglobal = atom_local_to_global[iatom][static_cast<std::size_t>(imu)];
                for (int jmu = 0; jmu != jnaux; ++jmu)
                {
                    const int jglobal = atom_local_to_global[jatom][static_cast<std::size_t>(jmu)];
                    block.payload.push_back(chi0_matrix[static_cast<std::size_t>(iglobal)
                                                            * static_cast<std::size_t>(naux)
                                                        + static_cast<std::size_t>(jglobal)]);
                }
            }
            blocks.push_back(std::move(block));
        }
    }

    const std::int32_t nblocks = checked_i32_from_size(blocks.size(), "Sternheimer chi0 v1 block count");
    std::int64_t offset = 6 * static_cast<std::int64_t>(sizeof(std::int32_t))
        + 2 * static_cast<std::int64_t>(sizeof(double))
        + static_cast<std::int64_t>(sizeof(std::int32_t))
        + static_cast<std::int64_t>(metadata.atom_naux.size() * sizeof(std::int32_t))
        + static_cast<std::int64_t>(blocks.size())
              * static_cast<std::int64_t>(sizeof(std::int32_t) + sizeof(std::int64_t));
    for (Chi0Block& block: blocks)
    {
        block.offset = offset;
        const std::size_t bytes = checked_mul_size(block.payload.size(),
                                                   sizeof(Complex),
                                                   "Sternheimer chi0 v1 block payload");
        offset += checked_i64_from_size(bytes, "Sternheimer chi0 v1 block payload");
    }

    static_assert(sizeof(Complex) == 2 * sizeof(double),
                  "Sternheimer chi0 v1 output expects complex<double> as two doubles.");
    std::ofstream out(filename.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
    if (!out.good())
    {
        throw std::runtime_error("Failed to open " + filename);
    }

    const std::int32_t iq = metadata.iq;
    const std::int32_t ifrequency = metadata.ifrequency;
    const std::int32_t naux_i32 = naux;
    const std::int32_t natom_i32 = checked_i32_from_size(natom, "Sternheimer chi0 v1 atom count");
    write_scalar(out, kChi0V1Marker, filename);
    write_scalar(out, iq, filename);
    write_scalar(out, ifrequency, filename);
    write_scalar(out, naux_i32, filename);
    write_scalar(out, kChi0V1ComplexFlag, filename);
    write_scalar(out, natom_i32, filename);
    write_scalar(out, metadata.omega, filename);
    write_scalar(out, metadata.weight, filename);
    write_scalar(out, nblocks, filename);
    for (const int atom_aux: metadata.atom_naux)
    {
        const std::int32_t atom_aux_i32 = atom_aux;
        write_scalar(out, atom_aux_i32, filename);
    }
    for (const Chi0Block& block: blocks)
    {
        const std::int32_t pair_index = block.pair_index;
        write_scalar(out, pair_index, filename);
        write_scalar(out, block.offset, filename);
    }
    for (const Chi0Block& block: blocks)
    {
        checked_write(out, block.payload.data(), block.payload.size() * sizeof(Complex), filename);
    }
    out.close();
}

SternheimerRPA::Complex SternheimerRPA::local_grid_dot(const Vector& lhs, const Vector& rhs, const double grid_weight)
{
    assert_same_size(lhs, rhs, "SternheimerRPA::local_grid_dot");
    Complex value(0.0, 0.0);
    for (std::size_t ir = 0; ir != lhs.size(); ++ir)
    {
        value += std::conj(lhs[ir]) * rhs[ir];
    }
    return grid_weight * value;
}

void SternheimerRPA::project_out_subspace(const std::vector<Vector>& subspace,
                                          const std::function<Complex(const Vector&, const Vector&)>& dot,
                                          Vector& vec)
{
    if (!dot)
    {
        throw std::invalid_argument("SternheimerRPA::project_out_subspace requires a dot callback.");
    }
    for (const Vector& basis_vec: subspace)
    {
        assert_same_size(basis_vec, vec, "SternheimerRPA::project_out_subspace");
        const Complex norm = dot(basis_vec, basis_vec);
        if (std::abs(norm) == 0.0)
        {
            continue;
        }
        const Complex coeff = dot(basis_vec, vec) / norm;
        for (std::size_t i = 0; i != vec.size(); ++i)
        {
            vec[i] -= coeff * basis_vec[i];
        }
    }
}

void SternheimerRPA::apply_kinetic_preconditioner(const std::vector<double>& kinetic_energy,
                                                  const double eigenvalue,
                                                  const double omega,
                                                  const double eta,
                                                  const Vector& input,
                                                  Vector& output)
{
    if (kinetic_energy.size() != input.size())
    {
        throw std::invalid_argument("SternheimerRPA::apply_kinetic_preconditioner size mismatch.");
    }
    output.resize(input.size());
    for (std::size_t ig = 0; ig != input.size(); ++ig)
    {
        const Complex denominator(kinetic_energy[ig] - eigenvalue + eta, omega);
        output[ig] = input[ig] / denominator;
    }
}

} // namespace ModuleRI
