#include "ions_move_cg.h"

#include "ions_move_basic.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
using namespace Ions_Move_Basic;

double Ions_Move_CG::RELAX_CG_THR = -1.0; // default is 0.5
//=================== NOTES ========================
// in vasp, it's like
// contolled by POTIM.
// the approximate minimum of the total energy
// is calculated from a cubic (or quadratic)
// interpolation taking into account the change
// of the total energy and the changes of forces.
// ( 3 pieces of information )
// 1. trial step
// 2. corrector step using cubic or quadratic interpolation
// 3. Check the force and see if we need brent method interpolation.
// 4. a new trial step....
// Brent's method is a complicated but popular
// root-finding algorithm combining the bisection method,
// the secant method and inverse quadratic interpolation
//=================== NOTES ========================

Ions_Move_CG::Ions_Move_CG()
{
    this->pos0 = nullptr;
    this->grad0 = nullptr;
    this->cg_grad0 = nullptr;
    this->move0 = nullptr;
}

Ions_Move_CG::~Ions_Move_CG()
{
    delete[] pos0;
    delete[] grad0;
    delete[] cg_grad0;
    delete[] move0;
}

void Ions_Move_CG::allocate(void)
{
    ModuleBase::TITLE("Ions_Move_CG", "allocate");
    assert(dim > 0);
    delete[] pos0;
    delete[] grad0;
    delete[] cg_grad0;
    delete[] move0;
    this->pos0 = new double[dim];
    this->grad0 = new double[dim];
    this->cg_grad0 = new double[dim];
    this->move0 = new double[dim];
    ModuleBase::GlobalFunc::ZEROS(pos0, dim);
    ModuleBase::GlobalFunc::ZEROS(grad0, dim);
    ModuleBase::GlobalFunc::ZEROS(cg_grad0, dim);
    ModuleBase::GlobalFunc::ZEROS(move0, dim);
    this->e0 = 0.0;
}

void Ions_Move_CG::start(UnitCell &ucell, const ModuleBase::matrix &force, const double &etot_in)
{
    ModuleBase::TITLE("Ions_Move_CG", "start");
    assert(dim > 0);
    assert(pos0 != 0);
    assert(grad0 != 0);
    assert(cg_grad0 != 0);
    assert(move0 != 0);

    static bool
        sd; // sd , trial are two parameters, when sd=trial=true ,a new direction begins, when sd = false trial =true
    static bool trial;  // a cubic interpolation is used to make the third point ,when sa = trial = false , we use Brent
                        // to get the
    static int ncggrad; // minimum point in this direction.
    static double fa, fb, fc; // ncggrad is a parameter to control the cg method , every ten cg direction , we change
                              // the direction back to
    static double xa, xb, xc, xpt, steplength, fmax; // the steepest descent method
    static int nbrent;

    double *pos = new double[dim];
    double *grad = new double[dim];
    double *cg_gradn = new double[dim];
    double *move = new double[dim];
    double *cg_grad = new double[dim];
    double best_x;
    double fmin;

    int flag = 0;

    ModuleBase::GlobalFunc::ZEROS(pos, dim);
    ModuleBase::GlobalFunc::ZEROS(grad, dim);
    ModuleBase::GlobalFunc::ZEROS(move, dim);
    ModuleBase::GlobalFunc::ZEROS(cg_grad, dim);

CG_begin:

    if (Ions_Move_Basic::istep == 1)
    {
        steplength = Ions_Move_Basic::relax_bfgs_init; // read in the init trust radius
        // std::cout<<"Ions_Move_Basic::relax_bfgs_init = "<<Ions_Move_Basic::relax_bfgs_init<<std::endl;
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

    Ions_Move_Basic::setup_gradient(ucell, force, pos, grad);
    // use energy_in and istep to setup etot and etot_old.
    Ions_Move_Basic::setup_etot(etot_in, 0);
    // use gradient and etot and etot_old to check
    // if the result is converged.

    // std::cout<<"sd = "<<sd<<"  trial = "<<trial<<"  istep = "<<istep<<std::endl;
    if (flag == 0)
    {
        Ions_Move_Basic::check_converged(ucell, grad);
        // std::cout<<"Ions_Move_Basic::converged = "<<Ions_Move_Basic::converged<<std::endl;
    }

    if (Ions_Move_Basic::converged)
    {
        Ions_Move_Basic::terminate(ucell);
    }
    else
    {
        if (sd)
        {
            e0 = etot_in;
            setup_cg_grad(grad,
                          grad0,
                          cg_grad,
                          cg_grad0,
                          ncggrad,
                          flag); // we use the last direction ,the last grad and the grad now to get the direction now
            ncggrad++;

            double norm = 0.0;
            for (int i = 0; i < dim; ++i)
            {
                norm += pow(cg_grad[i], 2);
            }
            norm = sqrt(norm);

            if (norm != 0.0)
            {
                for (int i = 0; i < dim; ++i)
                {
                    cg_gradn[i] = cg_grad[i] / norm;
                }
            }

            setup_move(move0, cg_gradn, steplength); // move the atom position
            Ions_Move_Basic::move_atoms(ucell, move0, pos);

            for (int i = 0; i < dim; i++) // grad0 ,cg_grad0 are used to store the grad and cg_grad for the future using
            {
                grad0[i] = grad[i];
                cg_grad0[i] = cg_grad[i];
            }

            f_cal(move0, move0, dim, xb); // xb = trial steplength
            f_cal(move0, grad, dim, fa);  // fa is the projection force in this direction

            fmax = fa;
            sd = false;

            if (GlobalV::RELAX_METHOD == "cg_bfgs")
            {
                if (Ions_Move_Basic::largest_grad * ModuleBase::Ry_to_eV / 0.529177
                    < RELAX_CG_THR) // cg to bfgs  by pengfei 13-8-8
                {
                    GlobalV::RELAX_METHOD = "bfgs";
                }
                Ions_Move_Basic::best_xxx = steplength;
            }

            Ions_Move_Basic::relax_bfgs_init = xb;
        }
        else
        {
            if (trial)
            {
                double e1 = etot_in;
                f_cal(move0, grad, dim, fb);
                f_cal(move0, move0, dim, xb);

                if ((std::abs(fb) < std::abs((fa) / 10.0)))
                {
                    sd = true;
                    trial = true;
                    steplength = xb;
                    flag = 1;
                    goto CG_begin;
                }

                double norm = 0.0;
                for (int i = 0; i < dim; ++i)
                {
                    norm += pow(cg_grad0[i], 2);
                }
                norm = sqrt(norm);

                if (norm != 0.0)
                {
                    for (int i = 0; i < dim; ++i)
                    {
                        cg_gradn[i] = cg_grad0[i] / norm;
                    }
                }

                third_order(e0, e1, fa, fb, xb, best_x); // cubic interpolation

                if (best_x > 6 * xb || best_x < (-xb))
                {
                    best_x = 6 * xb;
                }

                setup_move(move, cg_gradn, best_x);
                Ions_Move_Basic::move_atoms(ucell, move, pos);

                trial = false;
                xa = 0;
                f_cal(move0, move, dim, xc);
                xc = xb + xc;
                xpt = xc;
                Ions_Move_Basic::relax_bfgs_init = xc;
            }
            else
            {
                double xtemp, ftemp;
                f_cal(move0, grad, dim, fc);

                fmin = std::abs(fc);
                nbrent++;
                // std::cout<<"nbrent = "<<nbrent<<std::endl;
                // std::cout<<"xa = "<<xa<<" xb = "<<xb<<" xc = "<<xc<<" fa = "<<fa<<" fb = "<<fb<<" fc =
                // "<<fc<<std::endl;

                if ((fmin < std::abs((fmax) / 10.0)) || (nbrent > 3))
                {
                    nbrent = 0;
                    sd = true;
                    trial = true;
                    steplength = xpt;
                    flag = 1;
                    goto CG_begin;
                }
                else
                {
                    Brent(fa, fb, fc, xa, xb, xc, best_x, xpt); // Brent method
                    // std::cout<<"xc = "<<xc<<std::endl;
                    if (xc < 0)
                    {
                        sd = true;
                        trial = true;
                        steplength = xb;
                        flag = 2;
                        goto CG_begin;
                    }

                    double norm = 0.0;
                    for (int i = 0; i < dim; ++i)
                    {
                        norm += pow(cg_grad0[i], 2);
                    }
                    norm = sqrt(norm);

                    if (norm != 0.0)
                    {
                        for (int i = 0; i < dim; ++i)
                        {
                            cg_gradn[i] = cg_grad0[i] / norm;
                        }
                    }

                    setup_move(move, cg_gradn, best_x);
                    Ions_Move_Basic::move_atoms(ucell, move, pos);

                    Ions_Move_Basic::relax_bfgs_init = xc;
                }
            }
        }
    }

    delete[] cg_grad;
    delete[] grad;
    delete[] cg_gradn;
    delete[] pos;
    delete[] move;

    return;
}

void Ions_Move_CG::setup_cg_grad(double *grad,
                                 const double *grad0,
                                 double *cg_grad,
                                 const double *cg_grad0,
                                 const int &ncggrad,
                                 int &flag)
{
    ModuleBase::TITLE("Ions_Move_CG", "setup_cg_grad");
    assert(Ions_Move_Basic::istep > 0);
    double gamma;
    double cg0_cg, cg0_cg0, cg0_g;

    if (ncggrad % 10000 == 0 || flag == 2)
    {
        for (int i = 0; i < dim; i++)
        {
            cg_grad[i] = grad[i];
        }
    }
    else
    {
        double gp_gp = 0.0;  // grad_p.grad_p
        double gg = 0.0;     // grad.grad
        double g_gp = 0.0;   // grad_p.grad
        double cgp_gp = 0.0; // cg_grad_p.grad_p
        double cgp_g = 0.0;  // cg_grad_p.grad
        for (int i = 0; i < dim; i++)
        {
            gp_gp += grad0[i] * grad0[i];
            gg += grad[i] * grad[i];
            g_gp += grad0[i] * grad[i];
            cgp_gp += cg_grad0[i] * grad0[i];
            cgp_g += cg_grad0[i] * grad[i];
        }

        assert(g_gp != 0.0);
        const double gamma1 = gg / gp_gp; // FR
        // const double gamma2 = -(gg - g_gp)/(cgp_g - cgp_gp);  //CW
        const double gamma2 = (gg - g_gp) / gp_gp; // PRP
        // const double gamma = gg/cgp_gp;                      //D
        // const double gamma = -gg/(cgp_g - cgp_gp);             //D-Y

        if (gamma1 < 0.5)
        {
            gamma = gamma1;
        }
        else
        {
            gamma = gamma2;
        }

        for (int i = 0; i < dim; i++)
        {
            // we can consider step as modified gradient.
            cg_grad[i] = grad[i] + gamma * cg_grad0[i];
        }
    }
    return;
}

void Ions_Move_CG::third_order(const double &e0,
                               const double &e1,
                               const double &fa,
                               const double &fb,
                               const double x,
                               double &best_x)
{
    double k3, k2, k1;
    double dmoveh, dmove1, dmove2, dmove, ecal1, ecal2;

    k3 = 3 * ((fb + fa) * x - 2 * (e1 - e0)) / (x * x * x);
    k2 = (fb - fa) / x - k3 * x;
    k1 = fa;

    dmoveh = x * fb / (fa - fb);
    dmove1 = -k2 * (1 - sqrt(1 - 4 * k1 * k3 / (k2 * k2))) / (2 * k3);
    dmove2 = -k2 * (1 + sqrt(1 - 4 * k1 * k3 / (k2 * k2))) / (2 * k3);

    if ((std::abs(k3 / k1) < 0.01) || ((k1 * k3 / (k2 * k2)) >= 0.25)) // this condition may be wrong
    {
        dmove = dmoveh;
    }
    else
    {
        dmove1 = -k2 * (1 - sqrt(1 - 4 * k1 * k3 / (k2 * k2))) / (2 * k3);
        dmove2 = -k2 * (1 + sqrt(1 - 4 * k1 * k3 / (k2 * k2))) / (2 * k3);
        ecal1 = k3 * dmove1 * dmove1 * dmove1 / 3 + k2 * dmove1 * dmove1 / 2 + k1 * dmove1;
        ecal2 = k3 * dmove2 * dmove2 * dmove2 / 3 + k2 * dmove2 * dmove2 / 2 + k1 * dmove2;
        if (ecal2 > ecal1)
            dmove = dmove1 - x;
        else
            dmove = dmove2 - x;

        if (k3 < 0)
            dmove = dmoveh;
    }

    best_x = dmove;
    return;
}

void Ions_Move_CG::Brent(double &fa,
                         double &fb,
                         double &fc,
                         double &xa,
                         double &xb,
                         double &xc,
                         double &best_x,
                         double &xpt)
{
    double dmove;
    double tmp;
    double k2, k1, k0;
    double xnew1, xnew2;
    double ecalnew1, ecalnew2;

    if ((fa * fb) > 0)
    {
        dmove = (xc * fa - xa * fc) / (fa - fc);
        if (dmove > 4 * xc)
        // if(dmove > 4 * xc || dmove < 0)
        {
            dmove = 4 * xc;
        }
        xb = xc;
        fb = fc;
    }
    else
    {
        k2 = -((fb - fc) / (xb - xc) - (fa - fc) / (xa - xc)) / (xa - xb);
        k1 = (fa - fc) / (xa - xc) - k2 * (xa + xc);
        k0 = fa - k1 * xa - k2 * xa * xa;
        xnew1 = (-k1 - sqrt(k1 * k1 - 4 * k2 * k0)) / (2 * k2);
        xnew2 = (-k1 + sqrt(k1 * k1 - 4 * k2 * k0)) / (2 * k2);

        if (xnew1 > xnew2)
        {
            tmp = xnew2;
            xnew2 = xnew1;
            xnew1 = tmp;
        }

        ecalnew1 = k2 * xnew1 * xnew1 * xnew1 / 3 + k1 * xnew1 * xnew1 / 2 + k0 * xnew1;
        ecalnew2 = k2 * xnew2 * xnew2 * xnew2 / 3 + k1 * xnew2 * xnew2 / 2 + k0 * xnew2;
        dmove = xnew1;

        if (ecalnew1 > ecalnew2)
        {
            dmove = xnew2;
        }
        if (dmove < 0)
        {
            dmove = 2 * xc; // pengfei 14-6-5
        }
        if (fa * fc > 0)
        {
            xa = xc;
            fa = fc;
        }
        if (fb * fc > 0)
        {
            xb = xc;
            fb = fc;
        }
    }

    best_x = dmove - xpt;
    xpt = dmove;
    xc = dmove;

    return;
}

void Ions_Move_CG::f_cal(const double *g0, const double *g1, const int &dim, double &f_value)
{
    double hv0, hel;
    hel = 0;
    hv0 = 0;
    for (int i = 0; i < dim; i++)
    {
        hel += g0[i] * g1[i];
    }
    for (int i = 0; i < dim; i++)
    {
        hv0 += g0[i] * g0[i];
    }

    f_value = hel / sqrt(hv0);
    return;
}

void Ions_Move_CG::setup_move(double *move, double *cg_gradn, const double &trust_radius)
{
    // movement using gradient and trust_radius.
    for (int i = 0; i < dim; ++i)
    {
        move[i] = -cg_gradn[i] * trust_radius;
    }
    return;
}
