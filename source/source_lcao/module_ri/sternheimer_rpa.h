#ifndef STERNHEIMER_RPA_H
#define STERNHEIMER_RPA_H

#include <complex>
#include <functional>
#include <vector>

namespace ModuleRI
{

class SternheimerRPA
{
  public:
    using Complex = std::complex<double>;
    using Vector = std::vector<Complex>;

    struct SolverOptions
    {
        int max_iter = 300;
        double residual_tol = 1.0e-10;
        double breakdown_tol = 1.0e-14;
    };

    struct SolverResult
    {
        bool converged = false;
        int iterations = 0;
        double absolute_residual = 0.0;
        double relative_residual = 0.0;
    };

    struct LinearProblem
    {
        // Applies Pc (H - eps_i + i omega) Pc to a band-limited grid/PW vector.
        std::function<void(const Vector&, Vector&)> apply;

        // Applies a right preconditioner. If unset, identity is used.
        std::function<void(const Vector&, Vector&)> precondition;

        // Global inner product, including MPI reduction when running parallel.
        std::function<Complex(const Vector&, const Vector&)> dot;
    };

    static SolverResult solve_bicgstab(const LinearProblem& problem, const Vector& rhs, Vector& solution);

    static SolverResult solve_bicgstab(const LinearProblem& problem,
                                       const Vector& rhs,
                                       Vector& solution,
                                       const SolverOptions& options);

    static SolverResult solve_gmres(const LinearProblem& problem,
                                    const Vector& rhs,
                                    Vector& solution,
                                    const SolverOptions& options,
                                    int restart_dimension = 50);

    static void build_rhs_from_hartree_perturbation(const std::vector<double>& hartree_potential_r,
                                                    const Vector& psi_r,
                                                    Vector& rhs_r);

    static Complex accumulate_polarizability_grid_element(const std::vector<double>& hartree_potential_r,
                                                          const Vector& psi_r,
                                                          const Vector& delta_psi_r,
                                                          double grid_weight);

    static Complex local_grid_dot(const Vector& lhs, const Vector& rhs, double grid_weight);

    static void project_out_subspace(const std::vector<Vector>& subspace,
                                     const std::function<Complex(const Vector&, const Vector&)>& dot,
                                     Vector& vec);

    static void apply_kinetic_preconditioner(const std::vector<double>& kinetic_energy,
                                             double eigenvalue,
                                             double omega,
                                             double eta,
                                             const Vector& input,
                                             Vector& output);
};

} // namespace ModuleRI

#endif
