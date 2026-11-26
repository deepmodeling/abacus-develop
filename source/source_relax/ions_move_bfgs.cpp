#include "ions_move_bfgs.h"

#include "ions_move_basic.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"

#include <algorithm>
#include <cmath>

//============= MAP OF BFGS ===========================
// (1) start() -> BFGS_Basic::check_converged()
// -> restart_bfgs() -> bfgs_routine() -> save_bfgs()
// (2) restart_bfgs -> restore trust radius or reset_hessian()
// (3) bfgs_routine -> new_step() or interpolation
//============= MAP OF BFGS ===========================

Ions_Move_BFGS::Ions_Move_BFGS()
{
    // Default values for BFGS
    init_done = false;
}

Ions_Move_BFGS::~Ions_Move_BFGS() {};

void Ions_Move_BFGS::allocate()
{
    ModuleBase::TITLE("Ions_Move_BFGS", "init");
    if (init_done)
    {
        return;
    }
    this->allocate_basic();

    // initialize data members
    // be set in save_bfgs() function.
    this->save_flag = false;
    this->init_done = true;
    return;
}

void Ions_Move_BFGS::reset()
{
    if (init_done)
    {
        std::fill(pos.begin(), pos.end(), 0.0);
        std::fill(pos_p.begin(), pos_p.end(), 0.0);
        std::fill(grad.begin(), grad.end(), 0.0);
        std::fill(grad_p.begin(), grad_p.end(), 0.0);
        std::fill(move.begin(), move.end(), 0.0);
        std::fill(move_p.begin(), move_p.end(), 0.0);

        this->reset_hessian();
    }
    this->save_flag = false;
    this->tr_min_hit = false;

    Ions_Move_Basic::trust_radius = 0.0;
    Ions_Move_Basic::trust_radius_old = 0.0;
}

bool Ions_Move_BFGS::start(
    UnitCell &ucell,
    const ModuleBase::matrix &force,
    const double &energy_in,
    const int istep,
    int &update_iter,
    std::ofstream &ofs,
    std::vector<double> &etot_info,
    const Relax_Criteria &criteria)
{
    ModuleBase::TITLE("Ions_Move_BFGS", "start");

    Ions_Move_Basic::setup_gradient(ucell, force, pos.data(), grad.data(), ofs);
    Ions_Move_Basic::setup_etot(energy_in, istep, etot_info);
    int checked_iter = update_iter;
    const bool converged = Ions_Move_Basic::check_converged(
        ucell, this->grad.data(), checked_iter, ofs, etot_info, criteria.force_thr, criteria.force_thr_ev,
        criteria.out_level, criteria.test_relax_method);

    if (converged && etot_info[0] <= etot_info[1])
    {
        update_iter = checked_iter;
        Ions_Move_Basic::terminate(converged, update_iter, ucell, istep, ofs);
        return true;
    }
    else
    {
        // setup_gradient scales by lat0 for the shared convergence check.
        // The BFGS coordinates are in Bohr, so their energy gradient is in Ry/Bohr.
        for (double &g : grad)
        {
            g /= ucell.lat0;
        }
        this->restart_bfgs(update_iter, ofs, criteria.test_relax_method);
        this->bfgs_routine(istep, update_iter, ofs, etot_info, criteria.out_level, criteria.test_relax_method);
        this->save_bfgs();

        Ions_Move_Basic::move_atoms(ucell, move.data(), pos.data(), ofs, criteria.test_relax_method);
        // pos includes the applied symmetry projection and stays unwrapped for this step.
        for (std::size_t i = 0; i < pos.size(); ++i)
        {
            move_p[i] = pos[i] - pos_p[i];
        }
        return false;
    }
}

void Ions_Move_BFGS::restart_bfgs(int &update_iter, std::ofstream &ofs, const int test_relax_method)
{
    ModuleBase::TITLE("Ions_Move_BFGS", "restart_bfgs");
    using namespace Ions_Move_Basic;

    if (save_flag)
    {
        trust_radius_old = trust_radius;
        if (test_relax_method)
        {
            ModuleBase::GlobalFunc::OUT(ofs, "trust_radius_old (bohr)", trust_radius_old);
        }
    }
    else
    {
        std::fill(pos_p.begin(), pos_p.end(), 0.0);
        std::fill(grad_p.begin(), grad_p.end(), 0.0);
        std::fill(move_p.begin(), move_p.end(), 0.0);
        update_iter = 0;
        reset_hessian();
        tr_min_hit = false;
    }
}

void Ions_Move_BFGS::bfgs_routine(
    const int istep,
    int &update_iter,
    std::ofstream &ofs,
    std::vector<double> &etot_info,
    const std::string &out_level,
    const int test_relax_method)
{
    ModuleBase::TITLE("Ions_Move_BFGS", "bfgs_routine");
    using namespace Ions_Move_Basic;

    const int dim = grad.size();

    // etot_info[0] = etot (current total energy)
    // etot_info[1] = etot_p (previous total energy)
    // ediff = etot_info[0] - etot_info[1] (computed on demand)

    if (etot_info[0] > etot_info[1])
    {
        double dE0s = 0.0;
        for (int i = 0; i < dim; i++)
        {
            dE0s += this->grad_p[i] * this->move_p[i];
        }

        const double taken = std::sqrt(dot_func(move_p.data(), move_p.data(), dim));
        const double den = etot_info[0] - etot_info[1] - dE0s;
        if (dE0s < 0.0 && den > 0.0)
        {
            trust_radius = -0.5 * dE0s * taken / den;

            if (test_relax_method)
            {
                ModuleBase::GlobalFunc::OUT(ofs, "dE0s", dE0s);
                ModuleBase::GlobalFunc::OUT(ofs, "den", den);
                ModuleBase::GlobalFunc::OUT(ofs, "interpolated trust radius", trust_radius);
            }
        }
        else
        {
            trust_radius = 0.5 * taken;
            ofs << " quadratic interpolation is impossible." << std::endl;
        }

        etot_info[0] = etot_info[1];
        // Recover the accepted origin in the current periodic image.
        for (int i = 0; i < dim; ++i)
        {
            pos[i] -= move_p[i];
        }
        grad = grad_p;

        if (trust_radius < relax_bfgs_rmin)
        {
            ofs << "trust_radius = " << trust_radius << std::endl;
            ofs << "relax_bfgs_rmin = " << relax_bfgs_rmin << std::endl;
            ofs << "relax_bfgs_rmax = " << relax_bfgs_rmax << std::endl;
            ofs << " trust_radius < relax_bfgs_rmin, reset bfgs history." << std::endl;

            if (tr_min_hit)
            {
                ModuleBase::WARNING_QUIT("move_ions", "trust radius is too small! Break down.");
            }

            this->reset_hessian();

            for (int i = 0; i < dim; i++)
            {
                this->move[i] = -grad[i];
            }

            trust_radius = std::min(relax_bfgs_rmin, relax_bfgs_rmax);
            if (taken > 0.0)
            {
                trust_radius = std::min(trust_radius, 0.5 * taken);
            }
            tr_min_hit = true;
        }
        else
        {
            move = move_p;
            tr_min_hit = false;
        }
    }
    else if (etot_info[0] <= etot_info[1])
    {
        this->new_step(update_iter, ofs, etot_info, test_relax_method);
    }

    if (out_level == "ie")
    {
        std::cout << " BFGS TRUST (Bohr)    : " << trust_radius << std::endl;
    }

    ModuleBase::GlobalFunc::OUT(ofs, "istep", istep);
    ModuleBase::GlobalFunc::OUT(ofs, "update iteration", update_iter);

    const double norm = std::sqrt(dot_func(move.data(), move.data(), dim));

    if (norm == 0.0)
    {
        ModuleBase::WARNING_QUIT("Ions_Move_BFGS", "BFGS: move-length unreasonably short");
    }
    if (norm > trust_radius)
    {
        const double scale = trust_radius / norm;
        for (double &component : move)
        {
            component *= scale;
        }
    }

    return;
}
