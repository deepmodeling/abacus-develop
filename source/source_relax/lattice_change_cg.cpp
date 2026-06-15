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
}

bool Lattice_Change_CG::start(UnitCell &ucell, const ModuleBase::matrix &stress_in, const double &etot_in, std::ofstream& ofs, std::vector<double>& etot_info)
{
    ModuleBase::TITLE("Lattice_Change_CG", "start");

    assert(lat0.size() == static_cast<size_t>(dim));
    assert(grad0.size() == static_cast<size_t>(dim));
    assert(cg_grad0.size() == static_cast<size_t>(dim));
    assert(move0.size() == static_cast<size_t>(dim));

    static bool sd = false;
    static bool trial = false;
    static int ncggrad = 0;
    static double fa = 0.0;
    static double fb = 0.0;
    static double fc = 0.0;
    static double xa = 0.0;
    static double xb = 0.0;
    static double xc = 0.0;
    static double xpt = 0.0;
    static double steplength = 0.0;
    static double fmax = 0.0;
    static int nbrent = 0;

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
            steplength = Lattice_Change_Basic::lattice_change_ini;
            sd = true;
            trial = true;
            ncggrad = 0;
            fa = 0.0;
            fb = 0.0;
            fc = 0.0;
            xa = 0.0;
            xb = 0.0;
            xc = 0.0;
            xpt = 0.0;
            fmax = 0.0;
            nbrent = 0;
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

        if (sd)
        {
            e0 = etot_in;
            CG_Base::setup_cg_grad(dim, grad.data(), grad0.data(), cg_grad.data(), cg_grad0.data(), ncggrad, flag);
            ncggrad++;

            CG_Base::normalize(dim, cg_gradn.data(), cg_grad.data());
            CG_Base::setup_move(dim, move0.data(), cg_gradn.data(), steplength);
            Lattice_Change_Basic::change_lattice(ucell, move0.data(), lat.data());

            for (int i = 0; i < dim; i++)
            {
                grad0[i] = grad[i];
                cg_grad0[i] = cg_grad[i];
            }

            CG_Base::f_cal(dim, move0.data(), move0.data(), xb);
            CG_Base::f_cal(dim, move0.data(), grad.data(), fa);

            fmax = fa;
            sd = false;

            Lattice_Change_Basic::lattice_change_ini = xb;
            return false;
        }

        if (trial)
        {
            double e1 = etot_in;
            CG_Base::f_cal(dim, move0.data(), grad.data(), fb);
            CG_Base::f_cal(dim, move0.data(), move0.data(), xb);

            if ((std::abs(fb) < std::abs((fa) / 10.0)))
            {
                sd = true;
                trial = true;
                steplength = xb;
                flag = 1;
                continue;
            }

            CG_Base::normalize(dim, cg_gradn.data(), cg_grad0.data());
            CG_Base::third_order(e0, e1, fa, fb, xb, best_x);

            if (best_x > 6 * xb || best_x < (-xb))
            {
                best_x = 6 * xb;
            }

            CG_Base::setup_move(dim, move.data(), cg_gradn.data(), best_x);
            Lattice_Change_Basic::change_lattice(ucell, move.data(), lat.data());

            trial = false;
            xa = 0;
            CG_Base::f_cal(dim, move0.data(), move.data(), xc);
            xc = xb + xc;
            xpt = xc;

            Lattice_Change_Basic::lattice_change_ini = xc;
            return false;
        }

        double xtemp = 0.0;
        double ftemp = 0.0;
        CG_Base::f_cal(dim, move0.data(), grad.data(), fc);

        fmin = std::abs(fc);
        nbrent++;

        if ((fmin < std::abs((fmax) / 10.0)) || (nbrent > 3))
        {
            nbrent = 0;
            sd = true;
            trial = true;
            steplength = xpt;
            flag = 1;
            continue;
        }

        CG_Base::Brent(fa, fb, fc, xa, xb, xc, best_x, xpt);
        if (xc < 0)
        {
            sd = true;
            trial = true;
            steplength = xb;
            flag = 2;
            continue;
        }

        CG_Base::normalize(dim, cg_gradn.data(), cg_grad0.data());
        CG_Base::setup_move(dim, move.data(), cg_gradn.data(), best_x);
        Lattice_Change_Basic::change_lattice(ucell, move.data(), lat.data());

        Lattice_Change_Basic::lattice_change_ini = xc;
        return false;
    }
}