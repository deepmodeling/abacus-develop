#include "ions_move_bfgs.h"

#include <algorithm>
#include "source_io/module_parameter/parameter.h"
#include "ions_move_basic.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"

//============= MAP OF BFGS ===========================
// (1) start() -> BFGS_Basic::check_converged()
// -> restart_bfgs() -> bfgs_routine() -> save_bfgs()
// (2) restart_bfgs -> check_move() -> reset_hessian()
// (3) bfgs_routine -> new_step() or interpolation
//============= MAP OF BFGS ===========================

Ions_Move_BFGS::Ions_Move_BFGS()
{
    // Default values for BFGS
    init_done = false;
}

Ions_Move_BFGS::~Ions_Move_BFGS(){};

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

bool Ions_Move_BFGS::start(UnitCell& ucell, const ModuleBase::matrix& force, const double& energy_in, const int istep, int& update_iter, std::ofstream& ofs, std::vector<double>& etot_info)
{
    ModuleBase::TITLE("Ions_Move_BFGS", "start");

    if (first_step)
    {
        Ions_Move_Basic::setup_gradient(ucell, force, this->pos.data(), this->grad.data(), ofs);
        first_step = false;
    }
    else
    {
        std::vector<double> pos_tmp(3 * ucell.nat);
        Ions_Move_Basic::setup_gradient(ucell, force, pos_tmp.data(), this->grad.data(), ofs);
    }
    Ions_Move_Basic::setup_etot(energy_in, istep, etot_info);
    bool converged = Ions_Move_Basic::check_converged(ucell, this->grad.data(), update_iter, ofs, etot_info);

    if (converged)
    {
        Ions_Move_Basic::terminate(converged, update_iter, ucell, istep, ofs);
        return true;
    }
    else
    {
        this->restart_bfgs(ucell.lat0, update_iter, ofs);
        this->bfgs_routine(ucell.lat0, istep, update_iter, ofs, etot_info);
        this->save_bfgs();

        Ions_Move_Basic::move_atoms(ucell, move.data(), pos.data(), ofs);
        return false;
    }
}

void Ions_Move_BFGS::restart_bfgs(const double& lat0, int& update_iter, std::ofstream& ofs)
{
    ModuleBase::TITLE("Ions_Move_BFGS", "restart_bfgs");

    using namespace Ions_Move_Basic;

    const int dim = Ions_Move_Basic::dim;

    if (this->save_flag)
    {
        // (1) calculate the old trust radius
        trust_radius_old = 0.0;
        for (int i = 0; i < dim; i++)
        {
            // be careful! now the pos is *lat0 (Bohr)!!
            // bug(periodic boundary) trust_radius_old += (pos[i] - pos_p[i])*(pos[i] - pos_p[i]);
            trust_radius_old += this->move_p[i] * this->move_p[i];
        }
        trust_radius_old = sqrt(trust_radius_old);

        if (PARAM.inp.test_relax_method)
        {
            ModuleBase::GlobalFunc::OUT(ofs, "trust_radius_old (bohr)", trust_radius_old);
        }

        // (2)
        // normalize previous move, used in the case
        // calculate the previous movement of atoms. why I don't save it ??????????
        for (int i = 0; i < dim; i++)
        {
            // mohan add 2010-07-26.
            // there must be one of the two has the correct sign and value.
            this->move_p[i] = this->check_move(lat0, pos[i], pos_p[i]) / trust_radius_old;
            // std::cout << " " << std::setw(20) << move_p[i] << std::setw(20) << dpmin << std::endl;
        }
    }
    else
    {
        //    bfgs initialization
        std::fill(pos_p.begin(), pos_p.end(), 0.0);
            std::fill(grad_p.begin(), grad_p.end(), 0.0);
            std::fill(move_p.begin(), move_p.end(), 0.0);

        update_iter = 0;

        // set the trust radius old as the initial trust radius.
        trust_radius_old = relax_bfgs_init;
        this->reset_hessian();

        /*
        std::ifstream hess_file("hess_in");
        if(hess_file)
        {
            int rank1 = 0, rank2 = 0;
            hess_file >> rank1 >> rank2;
            if(rank1 == dim && rank2 == dim)
            {
                GlobalV::ofs_running << "\n Reading the approximate inverse hessian from file"<<std::endl;

                for(int i=0;i<dim;i++)
                {
                    for(int j=0;j<dim;j++)
                    {
                        hess_file >> inv_hess(i, j);
                    }
                }
            }
        }
        hess_file.close();
        */

        this->tr_min_hit = false;
    }
    return;
}

void Ions_Move_BFGS::bfgs_routine(const double& lat0, const int istep, int& update_iter, std::ofstream& ofs, std::vector<double>& etot_info)
{
    ModuleBase::TITLE("Ions_Move_BFGS", "bfgs_routine");
    using namespace Ions_Move_Basic;

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

        double den = etot_info[0] - etot_info[1] - dE0s;

        if (den > 1.0e-16)
        {
            trust_radius = -0.5 * dE0s * trust_radius_old / den;

            if (PARAM.inp.test_relax_method)
            {
                ModuleBase::GlobalFunc::OUT(ofs, "dE0s", dE0s);
                ModuleBase::GlobalFunc::OUT(ofs, "den", den);
                ModuleBase::GlobalFunc::OUT(ofs, "interpolated trust radius", trust_radius);
            }
        }
        else if (den <= 1.0e-16)
        {
            trust_radius = 0.5 * trust_radius_old;
            ofs << " quadratic interpolation is impossible." << std::endl;
        }

        etot_info[0] = etot_info[1];
        for (int i = 0; i < dim; i++)
        {
            this->pos[i] = pos_p[i];
            this->grad[i] = grad_p[i];
        }

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

            trust_radius = relax_bfgs_rmin;
            tr_min_hit = true;
        }
        else if (trust_radius >= relax_bfgs_rmin)
        {
            for (int i = 0; i < dim; i++)
            {
                this->move[i] = this->move_p[i] / trust_radius_old;
            }
            tr_min_hit = false;
        }
    }
    else if (etot_info[0] <= etot_info[1])
    {
        this->new_step(lat0, update_iter, ofs, etot_info);
    }

    if (PARAM.inp.out_level == "ie")
    {
        std::cout << " BFGS TRUST (Bohr)    : " << trust_radius << std::endl;
    }

    ModuleBase::GlobalFunc::OUT(ofs, "istep", istep);
    ModuleBase::GlobalFunc::OUT(ofs, "update iteration", update_iter);

    // combine the direction and move length now
    double norm = dot_func(this->move.data(), this->move.data(), dim);
    norm = sqrt(norm);

    if (norm < 1.0e-16)
    {
        ModuleBase::WARNING_QUIT("Ions_Move_BFGS", "BFGS: move-length unreasonably short");
    }
    else
    {
        // new move using trust_radius is
        // move / |move| * trust_radius (Bohr)
        for (int i = 0; i < dim; i++)
        {
            move[i] *= Ions_Move_Basic::trust_radius / norm;
        }
    }

    return;
}
