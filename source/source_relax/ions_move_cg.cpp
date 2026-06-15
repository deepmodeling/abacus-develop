#include "ions_move_cg.h"
#include "ions_move_basic.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include <vector>

double Ions_Move_CG::RELAX_CG_THR = -1.0; // default is 0.5

Ions_Move_CG::Ions_Move_CG()
{
    this->e0 = 0.0;
}

void Ions_Move_CG::allocate(const int dim)
{
    ModuleBase::TITLE("Ions_Move_CG", "allocate");
    assert(dim > 0);

    this->pos0.resize(dim, 0.0);
    this->grad0.resize(dim, 0.0);
    this->cg_grad0.resize(dim, 0.0);
    this->move0.resize(dim, 0.0);
    this->e0 = 0.0;

    this->sd = false;
    this->trial = false;
    this->ncggrad = 0;
    this->nbrent = 0;
    this->fa = 0.0;
    this->fb = 0.0;
    this->fc = 0.0;
    this->xa = 0.0;
    this->xb = 0.0;
    this->xc = 0.0;
    this->xpt = 0.0;
    this->steplength = 0.0;
    this->fmax = 0.0;
}

bool Ions_Move_CG::start(UnitCell &ucell, const ModuleBase::matrix &force, const double &etot_in, const int istep, int& update_iter, std::ofstream& ofs, std::vector<double>& etot_info, std::vector<std::string>& relax_method)
{
    ModuleBase::TITLE("Ions_Move_CG", "start");
    assert(Ions_Move_Basic::dim > 0);

    const int dim = Ions_Move_Basic::dim;
    std::vector<double> pos(dim, 0.0);
    std::vector<double> grad(dim, 0.0);
    std::vector<double> cg_gradn(dim, 0.0);
    std::vector<double> move(dim, 0.0);
    std::vector<double> cg_grad(dim, 0.0);
    double best_x = 0.0;
    double fmin = 0.0;
    int flag = 0;

    while (true)
    {
        if (istep == 1)
        {
            this->steplength = Ions_Move_Basic::relax_bfgs_init;
            this->sd = true;
            this->trial = true;
            this->ncggrad = 0;
            this->fa = 0.0;
            this->fb = 0.0;
            this->fc = 0.0;
            this->xa = 0.0;
            this->xb = 0.0;
            this->xc = 0.0;
            this->xpt = 0.0;
            this->fmax = 0.0;
            this->nbrent = 0;
        }

        Ions_Move_Basic::setup_gradient(ucell, force, pos.data(), grad.data(), ofs);
        Ions_Move_Basic::setup_etot(etot_in, 0, istep, ofs, etot_info);

        bool converged = false;
        if (flag == 0)
        {
            converged = Ions_Move_Basic::check_converged(ucell, grad.data(), update_iter, ofs, etot_info);
        }
        if (converged)
        {
            Ions_Move_Basic::terminate(converged, update_iter, ucell, istep, ofs);
            return true;
        }

        if (this->sd)
        {
            e0 = etot_in;
            CG_Base::setup_cg_grad(dim, grad.data(), grad0.data(), cg_grad.data(), cg_grad0.data(), this->ncggrad, flag);
            this->ncggrad++;

            CG_Base::normalize(dim, cg_gradn.data(), cg_grad.data());
            CG_Base::setup_move(dim, move0.data(), cg_gradn.data(), this->steplength);
            Ions_Move_Basic::move_atoms(ucell, move0.data(), pos.data(), ofs);

            for (int i = 0; i < dim; i++)
            {
                grad0[i] = grad[i];
                cg_grad0[i] = cg_grad[i];
            }

            CG_Base::f_cal(dim, move0.data(), move0.data(), this->xb);
            CG_Base::f_cal(dim, move0.data(), grad.data(), this->fa);
            this->fmax = this->fa;
            this->sd = false;

            if (relax_method[0] == "cg_bfgs")
            {
                if (Ions_Move_Basic::largest_grad * ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A
                    < RELAX_CG_THR)
                {
                    relax_method[0] = "bfgs";
                    relax_method[1] = "1";
                }
                Ions_Move_Basic::best_xxx = this->steplength;
            }

            Ions_Move_Basic::relax_bfgs_init = this->xb;
            return false;
        }

        if (this->trial)
        {
            double e1 = etot_in;
            CG_Base::f_cal(dim, move0.data(), grad.data(), this->fb);
            CG_Base::f_cal(dim, move0.data(), move0.data(), this->xb);

            if ((std::abs(this->fb) < std::abs((this->fa) / 10.0)))
            {
                this->sd = true;
                this->trial = true;
                this->steplength = this->xb;
                flag = 1;
                continue;
            }

            CG_Base::normalize(dim, cg_gradn.data(), cg_grad0.data());
            CG_Base::third_order(e0, e1, this->fa, this->fb, this->xb, best_x);

            if (best_x > 6 * this->xb || best_x < (-this->xb))
            {
                best_x = 6 * this->xb;
            }

            CG_Base::setup_move(dim, move.data(), cg_gradn.data(), best_x);
            Ions_Move_Basic::move_atoms(ucell, move.data(), pos.data(), ofs);
            this->trial = false;
            this->xa = 0;
            CG_Base::f_cal(dim, move0.data(), move.data(), this->xc);
            this->xc = this->xb + this->xc;
            this->xpt = this->xc;
            Ions_Move_Basic::relax_bfgs_init = this->xc;
            return false;
        }

        double xtemp = 0.0;
        double ftemp = 0.0;
        CG_Base::f_cal(dim, move0.data(), grad.data(), this->fc);
        fmin = std::abs(this->fc);
        this->nbrent++;

        if ((fmin < std::abs((this->fmax) / 10.0)) || (this->nbrent > 3))
        {
            this->nbrent = 0;
            this->sd = true;
            this->trial = true;
            this->steplength = this->xpt;
            flag = 1;
            continue;
        }

        CG_Base::Brent(this->fa, this->fb, this->fc, this->xa, this->xb, this->xc, best_x, this->xpt);
        if (this->xc < 0)
        {
            this->sd = true;
            this->trial = true;
            this->steplength = this->xb;
            flag = 2;
            continue;
        }

        CG_Base::normalize(dim, cg_gradn.data(), cg_grad0.data());
        CG_Base::setup_move(dim, move.data(), cg_gradn.data(), best_x);
        Ions_Move_Basic::move_atoms(ucell, move.data(), pos.data(), ofs);
        Ions_Move_Basic::relax_bfgs_init = this->xc;
        return false;
    }
}