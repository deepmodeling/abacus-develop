#include "ions_move_bfgs.h"

#include "ions_move_basic.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"

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
        return;
    this->allocate_basic();

    // initialize data members
    // be set in save_bfgs() function.
    this->save_flag = false;
    this->init_done = true;
    return;
}

void Ions_Move_BFGS::start(UnitCell& ucell, const ModuleBase::matrix& force, const double& energy_in)
{
    ModuleBase::TITLE("Ions_Move_BFGS", "start");

    // istep must be set eariler.

    // use force to setup gradient.
    Ions_Move_Basic::setup_gradient(ucell, force, this->pos, this->grad);
    // use energy_in and istep to setup etot and etot_old.
    Ions_Move_Basic::setup_etot(energy_in, 0);
    // use gradient and etot and etot_old to check
    // if the result is converged.
    Ions_Move_Basic::check_converged(ucell, this->grad);

    if (Ions_Move_Basic::converged)
    {
        Ions_Move_Basic::terminate(ucell);
    }
    else
    {
        // [ if new step ]
        // reset trust_radius_old.
        // [ if run from previous saved info ]
        // the BFGS file is read from previous run.
        // and the move_p is renormalized.
        this->restart_bfgs(ucell.lat0);

        //[ if etot>etot_p ]
        // interpolation
        //[ if etot<etot_p ]
        // calculate the new step -> the new move using hessian
        // matrix, and set the new trust radius.
        // [compute the move at last]
        this->bfgs_routine(ucell.lat0);

        // get prepared for the next try.
        // even if the energy is higher, we save the information.
        this->save_bfgs();

        Ions_Move_Basic::move_atoms(ucell, move, pos);
    }
    return;
}

void Ions_Move_BFGS::restart_bfgs(const double& lat0)
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

        if (GlobalV::test_relax_method)
        {
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "trust_radius_old (bohr)", trust_radius_old);
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
        //	bfgs initialization
        ModuleBase::GlobalFunc::ZEROS(pos_p, dim);
        ModuleBase::GlobalFunc::ZEROS(grad_p, dim);
        ModuleBase::GlobalFunc::ZEROS(move_p, dim);

        Ions_Move_Basic::update_iter = 0;

        // set the trust radius old as the initial trust radius.
        trust_radius_old = relax_bfgs_init;
        this->reset_hessian();

        /*
        std::ifstream hess_file("hess_in");
        if(hess_file)
        {
            int rank1,rank2;
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

void Ions_Move_BFGS::bfgs_routine(const double& lat0)
{
    ModuleBase::TITLE("Ions_Move_BFGS", "bfgs_routine");
    using namespace Ions_Move_Basic;

    // the bfgs algorithm starts here
    if (etot > etot_p)
    {
        // the previous step is rejected, line search goes on
        // we believe that we are in a correct direction, what we should do
        // in this case is to find a better step length until it is accepted

        // s: trust_radius
        // s': trust_radius_old
        // assume: E(s) = a*s*s + b*s + c( we use E(0), dE(0), E(s') )
        // E(s') = etot

        //-------------------------------------------
        // E(0) = etot_p
        // ==> c = etot_p

        //-------------------------------------------
        // Let's see how to calculate b.
        // dE(s) = (2*a*s + b) * ds
        // ==> dE(0) = b * ds
        // b = dE(0)/ds = grad_p

        //-------------------------------------------
        // Let's see how to calculate s.
        // dE(s)/ds = 2a*s + b = 0;
        // ==> s = - 0.5 * b / a

        //-------------------------------------------
        // Let's see how to calculate a:
        // E(s') = a*s'*s' + b*s' + etot_p = etot
        // ==> a*s'*s' = etot-etot_p-b*s'

        //-------------------------------------------
        // Let's see how to calculate final trust_radius
        // s = -0.5 * b / a = -0.5 * (b*s'*s') / (etot-etot_p-b*s')
        // s_min = - 0.5 * ( dE(0)*s'*s' ) / ( E(s') - E(0) - dE(0)*s' )

        // just like how we estimate the step length in CG
        // dE0s : b*y = averaged_grad * trust_radius_old
        double dE0s = 0.0;
        for (int i = 0; i < dim; i++)
        {
            // because dE(s)/dR(move_p) = Force(grad)
            // so dE = dR * grad
            dE0s += this->grad_p[i] * this->move_p[i];
        }

        double den = etot - etot_p - dE0s;

        if (den > 1.0e-16)
        {
            // get optimized trust radius
            trust_radius = -0.5 * dE0s * trust_radius_old / den;

            if (GlobalV::test_relax_method)
            {
                ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "dE0s", dE0s);
                ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "den", den);
                ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "interpolated trust radius", trust_radius);
            }
            // std::cout << " Formula : " << etot << " * s^2 + " << dE0s << " * s + " << etot_p << std::endl;
            // std::cout << " Lowest point : " << trust_radius << std::endl;
        }
        else if (den <= 1.0e-16)
        {
            // no quadratic interpolation is possible
            // then do is again, but smaller raidus.
            trust_radius = 0.5 * trust_radius_old;

            GlobalV::ofs_running << " quadratic interpolation is impossible." << std::endl;
        }
        // values from the last succeseful bfgs step are restored
        etot = etot_p;
        for (int i = 0; i < dim; i++)
        {
            this->pos[i] = pos_p[i];
            this->grad[i] = grad_p[i];
        }

        if (trust_radius < relax_bfgs_rmin)
        {
            // we are trapped in this case..., so the algorithim must be restart
            // the history is reset
            // xiaohui add 2013-03-17
            GlobalV::ofs_running << "trust_radius = " << trust_radius << std::endl;
            GlobalV::ofs_running << "relax_bfgs_rmin = " << relax_bfgs_rmin << std::endl;
            GlobalV::ofs_running << "relax_bfgs_rmax = " << relax_bfgs_rmax << std::endl;
            // xiaohui add 2013-03-17
            GlobalV::ofs_running << " trust_radius < relax_bfgs_rmin, reset bfgs history." << std::endl;

            if (tr_min_hit)
            {
                // the history has already been reset at the previous step
                // something is going wrong
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
            // old bfgs direction ( normalized ) is recovered
            for (int i = 0; i < dim; i++)
            {
                this->move[i] = this->move_p[i] / trust_radius_old;
            }
            tr_min_hit = false;
        }
    }
    else if (etot <= etot_p)
    {
        this->new_step(lat0);
    }

    if (GlobalV::OUT_LEVEL == "ie")
    {
        std::cout << " BFGS TRUST (Bohr)    : " << trust_radius << std::endl;
    }

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "istep", Ions_Move_Basic::istep);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "update iteration", Ions_Move_Basic::update_iter);

    // combine the direction and move length now
    double norm = dot_func(this->move, this->move, dim);
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
