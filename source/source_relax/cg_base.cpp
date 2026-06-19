#include "cg_base.h"
#include <cmath>
#include "source_base/global_function.h"

void CG_Base::setup_cg_grad(const int dim, double *grad, const double *grad0,
                            double *cg_grad, const double *cg_grad0,
                            const int &ncggrad, int &flag)
{
    ModuleBase::TITLE("CG_Base", "setup_cg_grad");
    double gamma = 0.0;

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

        const double gamma1 = gg / gp_gp; // FR
        const double gamma2 = (gg - g_gp) / gp_gp; // PRP

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
            cg_grad[i] = grad[i] + gamma * cg_grad0[i];
        }
    }
}

void CG_Base::setup_move(const int dim, double *move, double *cg_gradn, const double &trust_radius)
{
    ModuleBase::TITLE("CG_Base", "setup_move");
    for (int i = 0; i < dim; ++i)
    {
        move[i] = -cg_gradn[i] * trust_radius;
    }
}

void CG_Base::Brent(double &fa, double &fb, double &fc,
                    double &xa, double &xb, double &xc,
                    double &best_x, double &xpt)
{
    ModuleBase::TITLE("CG_Base", "Brent");
    double dmove = 0.0;
    double tmp = 0.0;
    double k2 = 0.0;
    double k1 = 0.0;
    double k0 = 0.0;
    double xnew1 = 0.0;
    double xnew2 = 0.0;
    double ecalnew1 = 0.0;
    double ecalnew2 = 0.0;

    if ((fa * fb) > 0)
    {
        dmove = (xc * fa - xa * fc) / (fa - fc);
        if (dmove > 4 * xc)
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
            dmove = 2 * xc;
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
}

void CG_Base::f_cal(const int dim, const double *g0, const double *g1, double &f_value)
{
    ModuleBase::TITLE("CG_Base", "f_cal");
    double hv0 = 0.0;
    double hel = 0.0;
    for (int i = 0; i < dim; i++)
    {
        hel += g0[i] * g1[i];
        hv0 += g0[i] * g0[i];
    }
    f_value = hel / sqrt(hv0);
}

void CG_Base::third_order(const double &e0, const double &e1,
                          const double &fa, const double &fb,
                          const double x, double &best_x)
{
    ModuleBase::TITLE("CG_Base", "third_order");
    double k3 = 0.0;
    double k2 = 0.0;
    double k1 = 0.0;
    double dmoveh = 0.0;
    double dmove1 = 0.0;
    double dmove2 = 0.0;
    double dmove = 0.0;
    double ecal1 = 0.0;
    double ecal2 = 0.0;

    k3 = 3 * ((fb + fa) * x - 2 * (e1 - e0)) / (x * x * x);
    k2 = (fb - fa) / x - k3 * x;
    k1 = fa;

    dmoveh = x * fb / (fa - fb);
    dmove1 = -k2 * (1 - sqrt(1 - 4 * k1 * k3 / (k2 * k2))) / (2 * k3);
    dmove2 = -k2 * (1 + sqrt(1 - 4 * k1 * k3 / (k2 * k2))) / (2 * k3);

    if ((std::abs(k3 / k1) < 0.01) || ((k1 * k3 / (k2 * k2)) >= 0.25))
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
        {
            dmove = dmove1 - x;
        }
        else
        {
            dmove = dmove2 - x;
        }

        if (k3 < 0)
        {
            dmove = dmoveh;
        }
    }

    best_x = dmove;
}

void CG_Base::normalize(const int dim, double *cg_gradn, const double *cg_grad)
{
    ModuleBase::TITLE("CG_Base", "normalize");
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
}