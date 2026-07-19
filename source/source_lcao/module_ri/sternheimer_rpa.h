#ifndef STERNHEIMER_RPA_H
#define STERNHEIMER_RPA_H

#include <complex>
#include <cstdint>
#include <functional>
#include <string>
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

    struct AuxiliaryChannel
    {
        int channel_index = -1;
        int atom_index = -1;
        int atom_local_index = -1;
    };

    struct Chi0V1Metadata
    {
        int iq = 1;
        int ifrequency = 1;
        double omega = 0.0;
        double weight = 1.0;
        std::vector<int> atom_naux;
    };

    struct CoulombV1Matrix
    {
        int iq = 0;
        std::vector<int> atom_naux;
        std::vector<Complex> values;
    };

    struct TransitionEnergyWindow
    {
        double emin_ha = 0.0;
        double emax_ha = 0.0;
    };

    struct FrequencyGrid
    {
        std::vector<double> omega_ha;
        std::vector<double> weights_ha;
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

    static void build_rhs_from_hartree_perturbation(const Vector& hartree_potential_r,
                                                    const Vector& psi_r,
                                                    Vector& rhs_r);

    static Complex accumulate_polarizability_grid_element(const std::vector<double>& hartree_potential_r,
                                                          const Vector& psi_r,
                                                          const Vector& delta_psi_r,
                                                          double grid_weight);

    static Complex accumulate_polarizability_grid_element(const Vector& hartree_potential_r,
                                                          const Vector& psi_r,
                                                          const Vector& delta_psi_r,
                                                          double grid_weight);

    static void accumulate_chi0_branch_column(const std::vector<std::vector<double>>& hartree_potentials_r,
                                              const Vector& psi_r,
                                              const Vector& delta_psi_r,
                                              double grid_weight,
                                              double occupation,
                                              int column_index,
                                              std::vector<Complex>& branch_matrix);

    static void accumulate_chi0_branch_column(const std::vector<Vector>& hartree_potentials_r,
                                              const Vector& psi_r,
                                              const Vector& delta_psi_r,
                                              double grid_weight,
                                              double occupation,
                                              int column_index,
                                              std::vector<Complex>& branch_matrix);

    static std::vector<Complex> symmetrize_chi0_imaginary_frequency(const std::vector<Complex>& branch_matrix,
                                                                    int num_channels);

    static std::int32_t chi0_v1_marker();

    static TransitionEnergyWindow transition_energy_window_from_eigenvalues_ry(
        const std::vector<double>& eigenvalues_ry,
        const std::vector<double>& occupations,
        double occupation_tolerance = 1.0e-8);

    static FrequencyGrid generate_greenx_minimax_frequency_grid(int nfreq, double emin_ha, double emax_ha);

    static FrequencyGrid read_frequency_grid_file(const std::string& filename, int expected_size);

    static int frequency_owner_rank(int ifrequency_zero_based, int mpi_ranks, int rank_shift = 0);

    static void write_chi0_v1_file(const std::string& filename,
                                   const Chi0V1Metadata& metadata,
                                   const std::vector<AuxiliaryChannel>& channels,
                                   const std::vector<Complex>& chi0_matrix);

    static CoulombV1Matrix read_coulomb_v1_files(const std::vector<std::string>& filenames);

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
