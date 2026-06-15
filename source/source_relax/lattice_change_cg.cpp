#include "lattice_change_cg.h"

#include "lattice_change_basic.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include <vector>

using namespace Lattice_Change_Basic;

Lattice_Change_CG::Lattice_Change_CG()
{
    this->e0 = 0.0;
}

void Lattice_Change_CG::allocate(void)
{
    ModuleBase::TITLE("Lattice_Change_CG", "allocate");
    assert(dim > 0);

    this->lat0.resize(dim, 0.0);
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

bool Lattice_Change_CG::start(UnitCell &ucell, const ModuleBase::matrix &stress_in, const double &etot_in, std::ofstream& ofs, std::vector<double>& etot_info)
{
    ModuleBase::TITLE("Lattice_Change_CG", "start");

    assert(lat0.size() == static_cast<size_t>(dim));
    assert(grad0.size() == static_cast<size_t>(dim));
    assert(cg_grad0.size() == static_cast<size_t>(dim));
    assert(move0.size() == static_cast<size_t>(dim));

    std::vector<double> lat(dim, 0.0);
    std::vector<double> grad(dim, 0.0);
    std::vector<double> cg_gradn(dim, 0.0);
    std::vector<double> move(dim, 0.0);
    std::vector<double> cg_grad(dim, 0.0);
    double best_x = 0.0;
    double fmin = 0.0;
    int flag = 0;

    while (true)
    {
        if (Lattice_Change_Basic::stress_step == 1)
        {
            this->steplength = Lattice_Change_Basic::lattice_change_ini;
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

        ModuleBase::matrix stress(stress_in);
        Lattice_Change_Basic::setup_gradient(ucell, lat.data(), grad.data(), stress);
        Lattice_Change_Basic::setup_etot(etot_in, 0, etot_info);

        bool converged = false;
        if (flag == 0)
        {
            converged = Lattice_Change_Basic::check_converged(ucell, stress, grad.data(), ofs);
        }

        if (converged)
        {
            Lattice_Change_Basic::terminate(converged, ofs);
            return true;
        }

        if (this->sd)
        {
            e0 = etot_in;
            CG_Base::setup_cg_grad(dim, grad.data(), grad0.data(), cg_grad.data(), cg_grad0.data(), this->ncggrad, flag);
            this->ncggrad++;

            CG_Base::normalize(dim, cg_gradn.data(), cg_grad.data());
            CG_Base::setup_move(dim, move0.data(), cg_gradn.data(), this->steplength);
            Lattice_Change_Basic::change_lattice(ucell, move0.data(), lat.data());

            for (int i = 0; i < dim; i++)
            {
                grad0[i] = grad[i];
                cg_grad0[i] = cg_grad[i];
            }

            CG_Base::f_cal(dim, move0.data(), move0.data(), this->xb);
            CG_Base::f_cal(dim, move0.data(), grad.data(), this->fa);

            this->fmax = this->fa;
            this->sd = false;

            Lattice_Change_Basic::lattice_change_ini = this->xb;
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
            Lattice_Change_Basic::change_lattice(ucell, move.data(), lat.data());

            this->trial = false;
            this->xa = 0;
            CG_Base::f_cal(dim, move0.data(), move.data(), this->xc);
            this->xc = this->xb + this->xc;
            this->xpt = this->xc;

            Lattice_Change_Basic::lattice_change_ini = this->xc;
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
        Lattice_Change_Basic::change_lattice(ucell, move.data(), lat.data());

        Lattice_Change_Basic::lattice_change_ini = this->xc;
        return false;
    }
}