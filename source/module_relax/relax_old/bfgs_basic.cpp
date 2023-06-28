#include "bfgs_basic.h"

#include "ions_move_basic.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
using namespace Ions_Move_Basic;

double BFGS_Basic::relax_bfgs_w1 = -1.0; // default is 0.01
double BFGS_Basic::relax_bfgs_w2 = -1.0; // defalut is 0.05

BFGS_Basic::BFGS_Basic()
{
    pos = nullptr;
    pos_p = nullptr;
    grad = nullptr;
    grad_p = nullptr;
    move = nullptr;
    move_p = nullptr;

    bfgs_ndim = 1;
}

BFGS_Basic::~BFGS_Basic()
{
    delete[] pos;
    delete[] pos_p;
    delete[] grad;
    delete[] grad_p;
    delete[] move;
    delete[] move_p;
}

void BFGS_Basic::allocate_basic(void)
{
    assert(dim > 0);

    delete[] pos;
    delete[] pos_p;
    delete[] grad;
    delete[] grad_p;
    delete[] move;
    delete[] move_p;

    pos = new double[dim];
    pos_p = new double[dim];
    grad = new double[dim];
    grad_p = new double[dim];
    move = new double[dim];
    move_p = new double[dim];

    ModuleBase::GlobalFunc::ZEROS(pos, dim);
    ModuleBase::GlobalFunc::ZEROS(grad, dim);
    ModuleBase::GlobalFunc::ZEROS(pos_p, dim);
    ModuleBase::GlobalFunc::ZEROS(grad_p, dim);
    ModuleBase::GlobalFunc::ZEROS(move, dim);
    ModuleBase::GlobalFunc::ZEROS(move_p, dim);

    // init inverse Hessien matrix.
    inv_hess.create(dim, dim);

    return;
}

void BFGS_Basic::update_inverse_hessian(const double &lat0)
{
    //  ModuleBase::TITLE("Ions_Move_BFGS","update_inverse_hessian");
    assert(dim > 0);

    double *s = new double[dim];
    double *y = new double[dim];
    ModuleBase::GlobalFunc::ZEROS(s, dim);
    ModuleBase::GlobalFunc::ZEROS(y, dim);

    for (int i = 0; i < dim; i++)
    {
        //      s[i] = this->pos[i] - this->pos_p[i];
        //		mohan update 2010-07-27
        s[i] = this->check_move(lat0, pos[i], pos_p[i]);
        s[i] *= lat0;

        y[i] = this->grad[i] - this->grad_p[i];
    }

    double sdoty = 0.0;
    for (int i = 0; i < dim; i++)
    {
        sdoty += s[i] * y[i];
    }
    if (std::abs(sdoty) < 1.0e-16)
    {
        GlobalV::ofs_running << " WARINIG: unexpected behaviour in update_inverse_hessian" << std::endl;
        GlobalV::ofs_running << " Resetting bfgs history " << std::endl;
        this->reset_hessian();
        return;
    }

    double *Hs = new double[dim];
    double *Hy = new double[dim];
    double *yH = new double[dim];
    ModuleBase::GlobalFunc::ZEROS(Hs, dim);
    ModuleBase::GlobalFunc::ZEROS(Hy, dim);
    ModuleBase::GlobalFunc::ZEROS(yH, dim);

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
            this->inv_hess(i, j)
                += 1.0 / sdoty * ((1.0 + ydotHy / sdoty) * s[i] * s[j] - (s[i] * yH[j] + Hy[i] * s[j]));
        }
    }

    delete[] s;
    delete[] y;
    delete[] Hs;
    delete[] Hy;
    delete[] yH;
    return;
}

void BFGS_Basic::check_wolfe_conditions(void)
{
    double dot_p = dot_func(grad_p, move_p, dim);
    double dot = dot_func(grad, move_p, dim);

    // if the total energy falls rapidly, enlarge the trust radius.
    bool wolfe1 = (etot - etot_p) < this->relax_bfgs_w1 * dot_p;

    // if the force is still very large, enlarge the trust radius,
    // otherwise the dot should be very small, in this case,
    // enlarge trst radius is not good.
    bool wolfe2 = std::abs(dot) > -this->relax_bfgs_w2 * dot_p;

    if (GlobalV::test_relax_method)
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "etot - etot_p", etot - etot_p);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "relax_bfgs_w1 * dot_p", relax_bfgs_w1 * dot_p);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "dot", dot);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "relax_bfgs_w2 * dot_p", relax_bfgs_w2 * dot_p);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "relax_bfgs_w1", relax_bfgs_w1);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "relax_bfgs_w2", relax_bfgs_w2);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "wolfe1", wolfe1);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "wolfe2", wolfe2);
    }

    this->wolfe_flag = wolfe1 && wolfe2;

    /*
       for(int i=0; i<dim; i++)
       {
       std::cout << " grad_p[" << i << "]=" << grad_p[i] << std::endl;
       }

       for(int i=0; i<dim; i++)
       {
       std::cout << " move_p[" << i << "]=" << move_p[i] << std::endl;
       }
     */

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "etot - etot_p", etot - etot_p);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "relax_bfgs_w1 * dot_p", relax_bfgs_w1 * dot_p);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "wolfe1", wolfe1);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "wolfe2", wolfe2);
    //  ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"dot = ",dot);
    //  ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"dot_p = ",dot_p);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "wolfe condition satisfied", wolfe_flag);
    return;
}

void BFGS_Basic::reset_hessian(void)
{
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            if (i == j)
                inv_hess(i, j) = 1.0;
            else
                inv_hess(i, j) = 0.0;
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
        this->move_p[i] = this->move[i];
    }
    return;
}

// a new bfgs step is done
// we have already done well in the previous direction
// we should get a new direction in this case
void BFGS_Basic::new_step(const double &lat0)
{
    ModuleBase::TITLE("BFGS_Basic", "new_step");

    //--------------------------------------------------------------------
    ++Ions_Move_Basic::update_iter;
    if (Ions_Move_Basic::update_iter == 1)
    {
        // Ions_Move_Basic::update_iter == 1 in this case
        // we haven't succes before, but we also need to decide a direction
        // this is the case when BFGS first start
        // if the gradient is very small now,
        // we don't need large relax_bfgs_init,
        // we choose a smaller one.
        if (Ions_Move_Basic::largest_grad < 0.01)
        {
            relax_bfgs_init = std::min(0.2, relax_bfgs_init);
        }

        Ions_Move_Basic::best_xxx = std::fabs(Ions_Move_Basic::best_xxx);

        // std::cout << "best_xxx=" << " " << best_xxx <<std::endl;

        relax_bfgs_init
            = std::min(Ions_Move_Basic::best_xxx, relax_bfgs_init); // cg to bfgs initial trust_radius   13-8-10 pengfei
    }
    else if (Ions_Move_Basic::update_iter > 1)
    {
        this->check_wolfe_conditions();
        this->update_inverse_hessian(lat0);
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
        GlobalV::ofs_running << " check the norm of new move " << dot_func(move, move, dim) << " (Bohr)" << std::endl;
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
        GlobalV::ofs_running << " Uphill move : resetting bfgs history" << std::endl;
        for (int i = 0; i < dim; i++)
        {
            move[i] = -grad[i];
        }
        this->reset_hessian();
    }

    //--------------------------------------------------------------------
    // the step must done after hessian is multiplied to grad.
    // std::cout<<"update_iter="<<Ions_Move_Basic::update_iter<<std::endl;
    if (Ions_Move_Basic::update_iter == 1)
    {
        trust_radius = relax_bfgs_init;

        this->tr_min_hit = false;
    }
    else if (Ions_Move_Basic::update_iter > 1)
    {
        trust_radius = trust_radius_old;
        this->compute_trust_radius();
    }
    // std::cout<<"trust_radius ="<<" "<<trust_radius;
    return;
}

// trust radius is computed in this function
// trust radius determine the step length
void BFGS_Basic::compute_trust_radius(void)
{
    ModuleBase::TITLE("BFGS_Basic", "compute_trust_radius");

    // (1) judge 1
    double dot = dot_func(grad_p, move_p, dim);
    bool ltest = (etot - etot_p) < this->relax_bfgs_w1 * dot;

    // (2) judge 2
    // calculate the norm of move, which
    // is used to compare to trust_radius_old.
    double norm_move = dot_func(this->move, this->move, dim);
    norm_move = std::sqrt(norm_move);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "move(norm)", norm_move);

    ltest = ltest && (norm_move > trust_radius_old);

    // (3) decide a
    double a;
    if (ltest)
    {
        a = 1.5;
    }
    else
    {
        a = 1.1;
    }

    /*
    std::cout << " a=" << a << std::endl;
    std::cout << " norm_move=" << norm_move << std::endl;
    std::cout << " trust_radius=" << trust_radius << std::endl;
    std::cout << " trust_radius_old=" << trust_radius_old << std::endl;
    */

    if (this->wolfe_flag)
    {
        trust_radius = std::min(relax_bfgs_rmax, 2.0 * a * trust_radius_old);
    }
    else
    {
        // mohan fix bug 2011-03-13 2*a*trust_radius_old -> a*trust_radius_old
        trust_radius = std::min(relax_bfgs_rmax, a * trust_radius_old);
        trust_radius = std::min(trust_radius, norm_move);
    }

    if (GlobalV::test_relax_method)
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "wolfe_flag", wolfe_flag);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "trust_radius_old", trust_radius_old);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "2*a*trust_radius_old", 2.0 * a * trust_radius_old);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "norm_move", norm_move);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Trust_radius (Bohr)", trust_radius);
    }

    if (trust_radius < relax_bfgs_rmin)
    {
        // the history should be reset, we got trapped
        if (tr_min_hit)
        {
            // the history has already been reset at the previous step
            // something is going wrongsomething is going wrong
            ModuleBase::WARNING_QUIT("bfgs", "bfgs history already reset at previous step, we got trapped!");
        }
        GlobalV::ofs_running << " Resetting BFGS history." << std::endl;
        this->reset_hessian();
        for (int i = 0; i < dim; i++)
        {
            move[i] = -grad[i];
        }
        trust_radius = relax_bfgs_rmin;
        tr_min_hit = true;
    }
    else
    {
        tr_min_hit = false;
    }

    return;
}

double BFGS_Basic::check_move(const double &lat0, const double &pos, const double &pos_p)
{
    // this must be careful.
    // unit is ucell.lat0.
    assert(lat0 > 0.0);
    const double direct_move = (pos - pos_p) / lat0;
    double shortest_move = direct_move;
    for (int cell = -1; cell <= 1; ++cell)
    {
        const double now_move = direct_move + cell;
        if (std::abs(now_move) < std::abs(shortest_move))
        {
            shortest_move = now_move;
        }
    }
    return shortest_move;
}
