// Wenfei Li, November 2022
// A new implementation of CG relaxation
#ifndef RELAX1_H
#define RELAX1_H

#include "line_search.h"
#include "source_base/matrix.h"
#include "source_base/matrix3.h"
#include "source_cell/unitcell.h"
#include "source_io/module_parameter/input_parameter.h"
#include <fstream>
#include <vector>

class Relax
{
  public:
    Relax() {};
    ~Relax() {};

    // prepare for relaxation
    void init_relax(const int nat_in, const Input_para& inp);
    // perform a single relaxation step
    bool relax_step(UnitCell& ucell,
                    const ModuleBase::matrix& force,
                    const ModuleBase::matrix& stress,
                    const double etot_in,
                    std::ofstream& ofs_running);

  private:
    int istep = 0; // count ionic step

    // setup gradient based on force and stress
    // constraints are considered here
    // also check if relaxation has converged
    // based on threshold in force & stress
    bool setup_gradient(const UnitCell& ucell, const ModuleBase::matrix& force, const ModuleBase::matrix& stress, std::ofstream& ofs_running);

    // Compute the ionic gradient from the force (in eV/Angstrom), honoring
    // per-atom move flags, and return the largest component. Records the value
    // in the force history and the running log, prints the optional out_level
    // diagnostics, and updates force_converged.
    double setup_ion_gradient(const UnitCell& ucell, const ModuleBase::matrix& force,
                              bool& force_converged, std::ofstream& ofs_running);

    // Compute the cell gradient from the stress, applying the fixed_axes
    // constraints (shape / volume / per-axis), and update force_converged.
    // Records the largest stress in the history and the running log. Only runs
    // when if_cell_moves is true.
    void setup_cell_gradient(const UnitCell& ucell, const ModuleBase::matrix& stress,
                             bool& force_converged, std::ofstream& ofs_running);

    // Print the converged / not-converged summary to the running log. On
    // convergence, reports the method, the converged step (istep + 1, since
    // setup_gradient runs before istep is incremented), and the per-step
    // largest force (and stress when the cell moves).
    void print_gradient_summary(bool force_converged, std::ofstream& ofs_running);

    // check whether previous line search is done
    bool check_line_search();

    // if line search not done : perform line search
    void perform_line_search(std::ofstream& ofs_running);

    // if line search done: find new search direction and make a trial move
    void new_direction(std::ofstream& ofs_running);

    // move ions and lattice vectors
    void move_cell_ions(UnitCell& ucell, const bool is_new_dir, std::ofstream& ofs_running);

    // Step 1 of move_cell_ions: update latvec along the cell search
    // direction, honoring per-axis free flags, the volume constraint and
    // fixed_ibrav. Saves latvec at the start of each CG step. Only runs when
    // if_cell_moves is true.
    void update_lattice(UnitCell& ucell, double fac, bool is_new_dir);

    // Steps 2 & 3 of move_cell_ions: compute the ionic displacement along the
    // ion search direction (Cartesian -> direct via the OLD GT), apply the
    // per-atom move flags and symmetry, then update taud/tau and print the
    // structure.
    void update_ion_positions(UnitCell& ucell, double fac, std::ofstream& ofs_running);

    // Steps 4 & 6 of move_cell_ions: refresh a1/a2/a3, omega and the
    // reciprocal lattice (G/GT/GGT) from the new latvec, broadcast them under
    // MPI, and re-setup the cell for the next SCF. Only runs when
    // if_cell_moves is true.
    void update_reciprocal_cell(UnitCell& ucell, std::ofstream& ofs_running);

    int nat = 0;         // number of atoms
    bool ltrial = false; // if last step is trial step

    double step_size = 0.0;

    // Gradients; _p means previous step
    ModuleBase::matrix grad_ion;
    ModuleBase::matrix grad_cell;
    ModuleBase::matrix grad_ion_p;
    ModuleBase::matrix grad_cell_p;

    // Search directions; _p means previous step
    ModuleBase::matrix search_dr_ion;
    ModuleBase::matrix search_dr_cell;
    ModuleBase::matrix search_dr_ion_p;
    ModuleBase::matrix search_dr_cell_p;

    // Used for applyting constraints
    bool if_cell_moves = false;

    // Keeps track of how many CG trial steps have been performed,
    // namely the number of CG directions followed
    // Note : this should not be confused with number of ionic steps
    // which includes both trial and line search steps
    int cg_step = 0;

    // in CG, search_dr = search_dr_p + grad * gamma
    double gamma = 0.0;
    void calculate_gamma();

    // Intermediate variables
    // I put them here because they are used across different subroutines
    double sr_sr = 0.0;
    double srp_srp = 0.0; // inner/cross products between search directions
    double gr_gr = 0.0;
    double gr_grp = 0.0;
    double grp_grp = 0.0; // inner/cross products between gradients
    double gr_sr = 0.0;   // cross product between search direction and gradient
    double e1ord1 = 0.0;
    double e1ord2 = 0.0;
    double e2ord = 0.0;
    double e2ord2 = 0.0;
    double dmove = 0.0;
    double dmovel = 0.0;
    double dmoveh = 0.0;
    double etot = 0.0;
    double etot_p = 0.0;
    /// previous cell volume in Angstrom^3, used to print volume diff during cell-relax
    double omega_p = 0.0;
    double force_thr_eva = 0.0;

    bool brent_done = false; // if brent line search is finished

    double fac_force = 0.0;
    double fac_stress = 0.0;

    ModuleBase::Matrix3 latvec_save;
    Line_Search ls;
    const Input_para* inp_ = nullptr;

    /// Largest force of each step (eV/Angstrom), includes trial and line-search steps.
    std::vector<double> max_force_history_;
    /// Largest stress of each step (kbar), only filled in cell-relax.
    std::vector<double> max_stress_history_;

  public:
    /// Read-only observers of the per-step convergence history, for the final summary.
    const std::vector<double>& get_max_force_history() const { return max_force_history_; }
    const std::vector<double>& get_max_stress_history() const { return max_stress_history_; }
};

#endif