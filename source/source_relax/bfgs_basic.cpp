#include "bfgs_basic.h"

#include "ions_move_basic.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"

#include <algorithm>
#include <cmath>

using namespace Ions_Move_Basic;

double BFGS_Basic::relax_bfgs_w1 = -1.0; // default is 0.01
double BFGS_Basic::relax_bfgs_w2 = -1.0; // default is 0.5

BFGS_Basic::BFGS_Basic() { bfgs_ndim = 1; }

void BFGS_Basic::allocate_basic(void)
{
    assert(dim > 0);

    pos.resize(dim, 0.0);
    pos_p.resize(dim, 0.0);
    grad.resize(dim, 0.0);
    grad_p.resize(dim, 0.0);
    move.resize(dim, 0.0);
    move_p.resize(dim, 0.0);

    // init inverse Hessien matrix.
    inv_hess.create(dim, dim);

    return;
}

void BFGS_Basic::update_inverse_hessian(bool first_update, std::ofstream &ofs)
{
    //  ModuleBase::TITLE("Ions_Move_BFGS","update_inverse_hessian");
    assert(dim > 0);

    const std::vector<double> &s = move_p;
    std::vector<double> y(dim, 0.0);

    for (int i = 0; i < dim; i++)
    {
        y[i] = this->grad[i] - this->grad_p[i];
    }

    double sdoty = 0.0;
    for (int i = 0; i < dim; i++)
    {
        sdoty += s[i] * y[i];
    }
    if (first_update)
    {
        // Set the initial inverse-Hessian scale from the first accepted secant pair.
        const double yy = dot_func(y.data(), y.data(), dim);
        if (yy > 0.0)
        {
            const double scale = std::abs(sdoty) / yy;
            if (scale > 0.0 && std::isfinite(scale))
            {
                inv_hess *= scale;
            }
        }
    }
    if (sdoty < 1.0e-16)
    {
        ofs << " Skipping inverse BFGS update with insufficient positive curvature." << std::endl;
        return;
    }

    std::vector<double> Hs(dim, 0.0);
    std::vector<double> Hy(dim, 0.0);
    std::vector<double> yH(dim, 0.0);

    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            Hs[i] += this->inv_hess(i, j) * s[j];
            Hy[i] += this->inv_hess(i, j) * y[j];
            yH[i] += y[j] * this->inv_hess(j, i); // mohan modify 2009-09-07
        }
    }

    double ydotHy = 0.0;
    for (int i = 0; i < dim; i++)
    {
        ydotHy += y[i] * Hy[i];
    }

    // inv_hess update
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            this->inv_hess(i, j) +=
                1.0 / sdoty * ((1.0 + ydotHy / sdoty) * s[i] * s[j] - (s[i] * yH[j] + Hy[i] * s[j]));
        }
    }

    return;
}

void BFGS_Basic::check_wolfe_conditions(std::ofstream &ofs, std::vector<double> &etot_info)
{
    double dot_p = dot_func(grad_p.data(), move_p.data(), dim);
    double dot = dot_func(grad.data(), move_p.data(), dim);

    // etot_info[0] = etot (current total energy)
    // etot_info[1] = etot_p (previous total energy)
    // ediff = etot_info[0] - etot_info[1] (computed on demand)
    const double ediff = etot_info[0] - etot_info[1];

    const bool wolfe1 = ediff <= relax_bfgs_w1 * dot_p;

    const bool wolfe2 = dot_p < 0.0 && std::abs(dot) <= -relax_bfgs_w2 * dot_p;

    this->wolfe_flag = wolfe1 && wolfe2;

    ModuleBase::GlobalFunc::OUT(ofs, "etot - etot_p", ediff);
    ModuleBase::GlobalFunc::OUT(ofs, "relax_bfgs_w1 * dot_p", relax_bfgs_w1 * dot_p);
    ModuleBase::GlobalFunc::OUT(ofs, "dot", dot);
    ModuleBase::GlobalFunc::OUT(ofs, "-relax_bfgs_w2 * dot_p", -relax_bfgs_w2 * dot_p);
    ModuleBase::GlobalFunc::OUT(ofs, "relax_bfgs_w1", relax_bfgs_w1);
    ModuleBase::GlobalFunc::OUT(ofs, "relax_bfgs_w2", relax_bfgs_w2);
    ModuleBase::GlobalFunc::OUT(ofs, "wolfe1", wolfe1);
    ModuleBase::GlobalFunc::OUT(ofs, "wolfe2", wolfe2);
    ModuleBase::GlobalFunc::OUT(ofs, "etot - etot_p", ediff);
    ModuleBase::GlobalFunc::OUT(ofs, "relax_bfgs_w1 * dot_p", relax_bfgs_w1 * dot_p);
    ModuleBase::GlobalFunc::OUT(ofs, "wolfe1", wolfe1);
    ModuleBase::GlobalFunc::OUT(ofs, "wolfe2", wolfe2);
    ModuleBase::GlobalFunc::OUT(ofs, "wolfe condition satisfied", wolfe_flag);
    return;
}

void BFGS_Basic::reset_hessian(void)
{
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            if (i == j)
            {
                inv_hess(i, j) = 1.0;
            }
            else
            {
                inv_hess(i, j) = 0.0;
            }
        }
    }
    return;
}

void BFGS_Basic::save_bfgs(void)
{
    this->save_flag = true;
    for (int i = 0; i < dim; i++)
    {
        this->pos_p[i] = this->pos[i];
        this->grad_p[i] = this->grad[i];
    }
    return;
}

// a new bfgs step is done
// we have already done well in the previous direction
// we should get a new direction in this case
void BFGS_Basic::new_step(
    int &update_iter,
    std::ofstream &ofs,
    std::vector<double> &etot_info,
    const int test_relax_method)
{
    ModuleBase::TITLE("BFGS_Basic", "new_step");

    //--------------------------------------------------------------------
    const int dim = grad.size();
    ++update_iter;
    if (update_iter > 1)
    {
        this->check_wolfe_conditions(ofs, etot_info);
        const bool first_update = update_iter == 2;
        this->update_inverse_hessian(first_update, ofs);
    }

    //--------------------------------------------------------------------
    // ---------------------------------
    // calculate the new move !!!
    // this step is very important !!!
    // ---------------------------------
    if (bfgs_ndim == 1)
    {
        // out.printrm("inv_hess",inv_hess,1.0e-8);
        for (int i = 0; i < dim; i++)
        {
            double tmp = 0.0;
            for (int j = 0; j < dim; j++)
            {
                tmp += this->inv_hess(i, j) * this->grad[j];
            }
            // we have got a new direction and step length!
            this->move[i] = -tmp;

            // std::cout << " move after hess " << move[i] << std::endl;
        }

        const double norm = std::sqrt(dot_func(move.data(), move.data(), dim));
        ofs << " check the norm of new move " << norm << " (Bohr)" << std::endl;
    }
    else if (bfgs_ndim > 1)
    {
        ModuleBase::WARNING_QUIT("Ions_Move_BFGS", "bfgs_ndim > 1 not implemented yet");
    }

    //--------------------------------------------------------------------
    // check our new direction here
    double dot = 0;
    for (int i = 0; i < dim; i++)
    {
        dot += grad[i] * move[i];
    }

    if (dot > 0.0)
    {
        ofs << " Uphill move : resetting bfgs history" << std::endl;
        for (int i = 0; i < dim; i++)
        {
            move[i] = -grad[i];
        }
        this->reset_hessian();
    }

    //--------------------------------------------------------------------
    // the step must done after hessian is multiplied to grad.
    // std::cout<<"update_iter="<<update_iter<<std::endl;
    if (update_iter == 1)
    {
        trust_radius = std::min({relax_bfgs_init, relax_bfgs_rmax, std::abs(best_xxx)});
        if (largest_grad < 0.01)
        {
            trust_radius = std::min(trust_radius, 0.2);
        }
        tr_min_hit = false;
    }
    else if (update_iter > 1)
    {
        this->compute_trust_radius(ofs, etot_info, test_relax_method);
    }
    // std::cout<<"trust_radius ="<<" "<<trust_radius;
    return;
}

// Grow the radius only after sufficient decrease on a step that reached the boundary.
void BFGS_Basic::compute_trust_radius(std::ofstream &ofs, std::vector<double> &etot_info, const int test_relax_method)
{
    ModuleBase::TITLE("BFGS_Basic", "compute_trust_radius");
    const int dim = move_p.size();
    const double taken = std::sqrt(dot_func(move_p.data(), move_p.data(), dim));
    const double slope = dot_func(grad_p.data(), move_p.data(), dim);
    const double ediff = etot_info[0] - etot_info[1];
    trust_radius = trust_radius_old;
    if (taken >= (1.0 - 1.e-10) * trust_radius_old && slope < 0.0 && ediff <= relax_bfgs_w1 * slope)
    {
        const double growth = wolfe_flag ? 2.0 : 1.5;
        trust_radius = std::min(relax_bfgs_rmax, growth * trust_radius_old);
    }
    tr_min_hit = false;
    if (test_relax_method)
    {
        ModuleBase::GlobalFunc::OUT(ofs, "wolfe_flag", wolfe_flag);
        ModuleBase::GlobalFunc::OUT(ofs, "step taken (Bohr)", taken);
        ModuleBase::GlobalFunc::OUT(ofs, "Trust_radius (Bohr)", trust_radius);
    }
}
