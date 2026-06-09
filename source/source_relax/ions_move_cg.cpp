#include "ions_move_cg.h"
#include "ions_move_basic.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include <vector>

using namespace Ions_Move_Basic;

double Ions_Move_CG::RELAX_CG_THR = -1.0; // default is 0.5

Ions_Move_CG::Ions_Move_CG()
{
    this->e0 = 0.0;
}

void Ions_Move_CG::allocate(void)
{
    ModuleBase::TITLE("Ions_Move_CG", "allocate");
    assert(dim > 0);

    this->pos0.resize(dim, 0.0);
    this->grad0.resize(dim, 0.0);
    this->cg_grad0.resize(dim, 0.0);
    this->move0.resize(dim, 0.0);
    this->e0 = 0.0;
}

void Ions_Move_CG::start(UnitCell &ucell, const ModuleBase::matrix &force, const double &etot_in)
{
    ModuleBase::TITLE("Ions_Move_CG", "start");
    assert(dim > 0);

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
        if (Ions_Move_Basic::istep == 1)
        {
            steplength = Ions_Move_Basic::relax_bfgs_init;
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

        Ions_Move_Basic::setup_gradient(ucell, force, pos.data(), grad.data());
        Ions_Move_Basic::setup_etot(etot_in, 0);

        if (flag == 0)
        {
            Ions_Move_Basic::check_converged(ucell, grad.data());
        }
        if (Ions_Move_Basic::converged)
        {
            Ions_Move_Basic::terminate(ucell);
            break;
        }

        if (sd)
        {
            e0 = etot_in;
            CG_Base::setup_cg_grad(dim, grad.data(), grad0.data(), cg_grad.data(), cg_grad0.data(), ncggrad, flag);
            ncggrad++;

            CG_Base::normalize(dim, cg_gradn.data(), cg_grad.data());
            CG_Base::setup_move(dim, move0.data(), cg_gradn.data(), steplength);
            Ions_Move_Basic::move_atoms(ucell, move0.data(), pos.data());

            for (int i = 0; i < dim; i++)
            {
                grad0[i] = grad[i];
                cg_grad0[i] = cg_grad[i];
            }

            CG_Base::f_cal(dim, move0.data(), move0.data(), xb);
            CG_Base::f_cal(dim, move0.data(), grad.data(), fa);
            fmax = fa;
            sd = false;

            if (Ions_Move_Basic::relax_method[0] == "cg_bfgs")
            {
                if (Ions_Move_Basic::largest_grad * ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A
                    < RELAX_CG_THR)
                {
                    Ions_Move_Basic::relax_method[0] = "bfgs";
                    Ions_Move_Basic::relax_method[1] = "1";
                }
                Ions_Move_Basic::best_xxx = steplength;
            }

            Ions_Move_Basic::relax_bfgs_init = xb;
            break;
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
            Ions_Move_Basic::move_atoms(ucell, move.data(), pos.data());
            trial = false;
            xa = 0;
            CG_Base::f_cal(dim, move0.data(), move.data(), xc);
            xc = xb + xc;
            xpt = xc;
            Ions_Move_Basic::relax_bfgs_init = xc;
            break;
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
        Ions_Move_Basic::move_atoms(ucell, move.data(), pos.data());
        Ions_Move_Basic::relax_bfgs_init = xc;
        break;
    }

    return;
}