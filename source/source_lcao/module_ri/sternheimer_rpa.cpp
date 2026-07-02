#include "source_lcao/module_ri/sternheimer_rpa.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace ModuleRI
{

namespace
{

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

} // namespace

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
