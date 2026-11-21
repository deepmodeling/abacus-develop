#include "relax_sync.h"

#include "relax_geometry.h"
#include "relax_history.h"
#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_base/tool_title.h"
#include "source_cell/print_cell.h"
#include "source_io/module_parameter/input_parameter.h"

#include <cmath>
#include <iostream>
#include <utility>

void Relax::init_relax(const int nat_in, const Input_para &inp)
{
    ModuleBase::TITLE("Relax", "init_relax");

    inp_ = &inp;

    max_force_history_.clear();
    max_stress_history_.clear();

    // set some initial conditions / constants
    nat = nat_in;
    istep = 0;
    cg_step = 0;
    ltrial = false;
    brent_done = false;
    step_size = 1.0;
    srp_srp = 100000;
    etot = 0;
    etot_p = 0;
    omega_p = 0.0;

    force_thr_eva = inp_->force_thr * ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A; // convert to eV/Angstrom
    fac_force = inp_->relax_scale_force * 0.1;
    fac_stress = fac_force / nat;

    // allocate some data structures

    // gradients in ion and cell; current and previous
    grad_ion.create(nat, 3);
    grad_cell.create(3, 3);

    grad_ion_p.create(nat, 3);
    grad_cell_p.create(3, 3);

    // search directions for ion and cell; current and previous
    search_dr_ion.create(nat, 3);
    search_dr_cell.create(3, 3);

    search_dr_ion_p.create(nat, 3);
    search_dr_cell_p.create(3, 3);

    // set if we are allowing lattice vectors to move
    if_cell_moves = false;
    if (inp_->calculation == "cell-relax")
    {
        if_cell_moves = true;
    }
}

bool Relax::relax_step(
    UnitCell &ucell,
    const ModuleBase::matrix &force,
    const ModuleBase::matrix &stress,
    const double etot_in,
    std::ofstream &ofs_running)
{
    ModuleBase::TITLE("Relax", "relax_step");

    ucell.ionic_position_updated = false;
    ucell.cell_parameter_updated = false;

    etot = etot_in * ModuleBase::Ry_to_eV; // convert to eV
    if (istep == 0)
    {
        etot_p = etot;
        omega_p = ucell.omega * pow(ModuleBase::BOHR_TO_A, 3);
    }

    bool relax_done = this->setup_gradient(ucell, force, stress, ofs_running);
    if (relax_done)
    {
        return relax_done;
    }

    this->calculate_gamma();

    bool ls_done = this->check_line_search();

    if (ls_done)
    {
        this->new_direction(ofs_running);
        this->move_cell_ions(ucell, true, ofs_running);
    }
    else
    {
        this->perform_line_search(ofs_running);
        this->move_cell_ions(ucell, false, ofs_running);
        dmovel = dmove;
    }

    istep++;

    return relax_done;
}

bool Relax::setup_gradient(
    const UnitCell &ucell,
    const ModuleBase::matrix &force,
    const ModuleBase::matrix &stress,
    std::ofstream &ofs_running)
{
    ModuleBase::TITLE("Relax", "setup_gradient");

    // if not relax, then return converged
    if (!(inp_->calculation == "relax" || inp_->calculation == "cell-relax"))
    {
        return true;
    }

    // indicating whether force & stress are converged
    bool force_converged = true;

    //=========================================
    // set gradient for ions degrees of freedom
    //=========================================

    setup_ion_gradient(ucell, force, force_converged, ofs_running);

    //=========================================
    // set gradient for cell degrees of freedom
    //=========================================

    if (if_cell_moves)
    {
        setup_cell_gradient(ucell, stress, force_converged, ofs_running);
    }

    print_gradient_summary(force_converged, ofs_running);

    return force_converged;
}

void Relax::setup_ion_gradient(
    const UnitCell &ucell,
    const ModuleBase::matrix &force,
    bool &force_converged,
    std::ofstream &ofs_running)
{
    const double max_grad = relax_geometry::project_ionic_forces(ucell, force, grad_ion);
    if (max_grad > force_thr_eva)
    {
        force_converged = false;
    }
    if (inp_->out_level == "ie")
    {
        if (if_cell_moves)
        {
            const double omega_ang = ucell.omega * pow(ModuleBase::BOHR_TO_A, 3);
            const double omega_diff = omega_ang - omega_p;
            const double omega_ratio = (std::abs(omega_p) > 0.0) ? omega_diff / omega_p * 100.0 : 0.0;
            std::cout << " CELL VOLUME (Angstroms^3)   : " << omega_ang << std::endl;
            std::cout << " VOLUME DIFF (Angstroms^3)   : " << omega_diff << std::endl;
            std::cout << " VOLUME RATIO (%)            : " << omega_ratio << std::endl;
            omega_p = omega_ang;
        }
        std::cout << " ETOT DIFF (eV)              : " << etot - etot_p << std::endl;
        std::cout << " LARGEST GRAD (eV/Angstrom)  : " << max_grad << std::endl;
        etot_p = etot;
    }

    max_force_history_.push_back(max_grad);
    ofs_running << "\n Largest force is " << max_grad << " eV/Angstrom while threshold is " << inp_->force_thr_ev
                << " eV/Angstrom" << std::endl;
}

void Relax::setup_cell_gradient(
    const UnitCell &ucell,
    const ModuleBase::matrix &stress,
    bool &force_converged,
    std::ofstream &ofs_running)
{
    const double largest_grad = relax_geometry::project_cell_stress(ucell, stress, inp_->fixed_axes, grad_cell);
    if (largest_grad > inp_->stress_thr)
    {
        force_converged = false;
    }

    max_stress_history_.push_back(largest_grad);
    ofs_running << " Largest stress is " << largest_grad << " kbar while threshold is " << inp_->stress_thr << " kbar"
                << std::endl;
}

void Relax::print_gradient_summary(bool force_converged, std::ofstream &ofs_running)
{
    if (force_converged)
    {
        // setup_gradient runs before istep is incremented at the end of
        // relax_step, so the step that just converged is istep + 1.
        const int converged_step = istep + 1;
        const std::string method = inp_->relax_method.empty() ? "unknown" : inp_->relax_method;
        ofs_running << " Relaxation method: " << method << std::endl;
        ofs_running << " Relaxation converged in " << converged_step << " step(s)." << std::endl;
        ofs_running << " Largest force per step (eV/Angstrom):" << format_relax_history(max_force_history_);
        if (if_cell_moves)
        {
            ofs_running << " Largest stress per step (kbar):" << format_relax_history(max_stress_history_);
        }
        ofs_running << "\n Relaxation is converged!" << std::endl;
    }
    else
    {
        ofs_running << "\n Relaxation is not converged yet!" << std::endl;
    }
}

void Relax::calculate_gamma()
{
    ModuleBase::TITLE("Relax", "calculate_gamma");

    // no need to calculate gamma if last step is trial
    // since we won't update search direction
    if (ltrial)
    {
        return;
    }

    grp_grp = 0.0; // grad_p*grad_p
    gr_grp = 0.0;  // grad  *grad_p
    gr_gr = 0.0;   // grad  *grad
    gr_sr = 0.0;   // grad  *search_dir

    for (int iat = 0; iat < nat; iat++)
    {
        for (int i = 0; i < 3; i++)
        {
            grp_grp += grad_ion_p(iat, i) * grad_ion_p(iat, i);
            gr_grp += grad_ion_p(iat, i) * grad_ion(iat, i);
            gr_gr += grad_ion(iat, i) * grad_ion(iat, i);
            gr_sr += grad_ion(iat, i) * search_dr_ion(iat, i);
        }
    }

    if (if_cell_moves)
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                grp_grp += grad_cell_p(i, j) * grad_cell_p(i, j) / nat;
                gr_grp += grad_cell_p(i, j) * grad_cell(i, j) / nat;
                gr_gr += grad_cell(i, j) * grad_cell(i, j) / nat;
                gr_sr += grad_cell(i, j) * search_dr_cell(i, j) / nat;
            }
        }
    }

    if (cg_step == 0)
    {
        gamma = 0.0;
    }
    else
    {
        gamma = (gr_gr - gr_grp) / grp_grp; // Polak-Riebere
    }
}

bool Relax::check_line_search()
{
    ModuleBase::TITLE("Relax", "check_line_search");

    // if last step is trial step towards new direction
    // then line search is not finished
    // we will perform line search
    if (ltrial)
    {
        ltrial = false;
        return false;
    }

    if (std::abs(gr_sr) * std::max(gamma, 1.0) > std::abs(gr_gr) / 5.0 &&
        !brent_done) // last brent line search not finished
    {
        return false;
    }

    return true;
}

void Relax::perform_line_search(std::ofstream &ofs_running)
{
    ModuleBase::TITLE("Relax", "line_search");

    double f = 0.0; // 1st order energy difference
    for (int iat = 0; iat < nat; iat++)
    {
        for (int i = 0; i < 3; i++)
        {
            f -= step_size * fac_force * search_dr_ion(iat, i) * grad_ion(iat, i);
        }
    }

    if (if_cell_moves)
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                f -= step_size * fac_stress * search_dr_cell(i, j) * grad_cell(i, j);
            }
        }
    }

    // perform line search
    bool restart_brent = false;
    double x = dmovel, y = etot;
    double xnew = 0.0, yd = 0.0;

    brent_done = this->ls.line_search(restart_brent, x, y, f, xnew, force_thr_eva, ofs_running);
    dmove = xnew;

    return;
}

void Relax::new_direction(std::ofstream &ofs_running)
{
    ModuleBase::TITLE("Relax", "new_direction");
    if (cg_step != 0)
    {
        step_size += 0.2 * step_size * (dmovel - 1.0);
    }

    // set GAMMA to zero if line minimization was not sufficient
    if (5.0 * std::abs(gr_sr) * gamma > std::abs(gr_gr))
    {
        gamma = 0.0;
        cg_step = 0; // reset cg
    }

    // perform trial step

    // set search vectors
    search_dr_ion_p = search_dr_ion;
    search_dr_cell_p = search_dr_cell;

    grad_ion_p = grad_ion;
    grad_cell_p = grad_cell;

    search_dr_ion = grad_ion + gamma * search_dr_ion;
    search_dr_cell = grad_cell + gamma * search_dr_cell;

    // modify step if necessary
    sr_sr = 1.0e-10;
    for (int iat = 0; iat < nat; iat++)
    {
        for (int i = 0; i < 3; i++)
        {
            sr_sr += 1 / fac_force * search_dr_ion(iat, i) * search_dr_ion(iat, i);
        }
    }

    if (if_cell_moves)
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                sr_sr += 1 / fac_stress * search_dr_cell(i, j) * search_dr_cell(i, j);
            }
        }
    }

    // if length of search vector increased, rescale step to avoid too large trial steps
    if (sr_sr > srp_srp)
    {
        step_size *= srp_srp / sr_sr;
    }
    srp_srp = sr_sr;

    double f = 0.0; // first order change in energy (gradient in the search direction)
    for (int iat = 0; iat < nat; iat++)
    {
        for (int i = 0; i < 3; i++)
        {
            f -= step_size * fac_force * search_dr_ion(iat, i) * grad_ion(iat, i);
        }
    }

    if (if_cell_moves)
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                f -= step_size * fac_stress * search_dr_cell(i, j) * grad_cell(i, j);
            }
        }
    }

    // prepare for line search
    bool restart = true;
    double x = 0, y = etot;
    // TODO: add a certain threshold for the progress.
    double xnew, yd = 1e-8;

    this->ls.line_search(restart, x, y, f, xnew, yd, ofs_running);

    dmovel = 1.0;
    ltrial = true;
    cg_step++;

    return;
}

void Relax::move_cell_ions(UnitCell &ucell, const bool is_new_dir, std::ofstream &ofs_running)
{
    ModuleBase::TITLE("Relax", "move_cell_ions");

    ucell.ionic_position_updated = true;
    const double fraction = is_new_dir ? 1.0 : dmove - dmovel;

    // Fractional displacements use the old GT. Cartesian positions use the updated lattice.
    const ModuleBase::Matrix3 old_gt = ucell.GT;
    if (if_cell_moves)
    {
        ucell.cell_parameter_updated = true;
        const ModuleBase::Matrix3 lattice_displacement = lattice_step(ucell, fraction, is_new_dir);
        const bool fixed_volume = inp_->fixed_axes == "volume";
        relax_geometry::apply_lattice_step(ucell, lattice_displacement, fixed_volume, inp_->fixed_ibrav);
    }
    ModuleBase::matrix ionic_displacement = search_dr_ion * (step_size * fac_force);
    relax_geometry::apply_ionic_step(ucell, old_gt, std::move(ionic_displacement), fraction);
    unitcell::print_tau(ucell.atoms, ucell.Coordinate, ucell.ntype, ucell.lat0, ofs_running);

    if (if_cell_moves)
    {
        relax_geometry::refresh_cell(ucell, inp_->nspin, ofs_running);
    }
}

ModuleBase::Matrix3 Relax::lattice_step(const UnitCell &ucell, double fraction, bool is_new_dir)
{
    // CG also uses the symmetrized, transposed search tensor in subsequent dot products.
    const bool use_symmetry = ModuleSymmetry::Symmetry::symm_flag && ucell.symm.nrotk > 0;
    if (use_symmetry)
    {
        search_dr_cell = ModuleBase::transpose(search_dr_cell);
        ucell.symm.symmetrize_mat3(search_dr_cell, ucell.lat);
    }

    const ModuleBase::matrix &dir = search_dr_cell;
    ModuleBase::Matrix3 direction(
        dir(0, 0), dir(0, 1), dir(0, 2), dir(1, 0), dir(1, 1), dir(1, 2), dir(2, 0), dir(2, 1), dir(2, 2));
    if (use_symmetry)
    {
        direction = direction.Transpose();
    }

    // Every correction along a line uses the lattice at the start of that CG direction.
    if (is_new_dir)
    {
        latvec_save = ucell.latvec;
    }
    return (latvec_save * direction) * (step_size * fraction * fac_stress);
}
