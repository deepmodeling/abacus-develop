#include "ions_move_sd.h"

#include <algorithm>
#include "source_io/module_parameter/parameter.h"
#include "ions_move_basic.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"

using namespace Ions_Move_Basic;

Ions_Move_SD::Ions_Move_SD() : energy_saved(1.0e10)
{
}

void Ions_Move_SD::allocate()
{
    ModuleBase::TITLE("Ions_Move_SD", "allocate");
    assert(dim > 0);
    grad_saved.resize(dim, 0.0);
    pos_saved.resize(dim, 0.0);
}

void Ions_Move_SD::start(UnitCell& ucell, const ModuleBase::matrix& force, const double& etot_in, std::ofstream& ofs)
{
    ModuleBase::TITLE("Ions_Move_SD", "start");

    assert(dim > 0);
    assert(grad_saved.size() == static_cast<size_t>(dim));
    assert(pos_saved.size() == static_cast<size_t>(dim));

    std::vector<double> pos(dim, 0.0);
    std::vector<double> grad(dim, 0.0);
    std::vector<double> move(dim, 0.0);

    // 1: ediff = 0
    // 0: ediff < 0
    bool judgement = false;
    setup_etot(etot_in, judgement);
    setup_gradient(ucell, force, pos.data(), grad.data(), ofs);

    if (istep == 1 || etot_in <= energy_saved)
    {
        printf("in cheak_converged");
        printf("pos[0]: %f\n", pos[0]);
        energy_saved = etot_in;
        for (int i = 0; i < dim; i++) 
        {
            pos_saved[i] = pos[i];
        }
        for (int i = 0; i < dim; i++)
        {
            grad_saved[i] = grad[i];
        }
        // normalize the gradient, in convinience to
        // move atom.
        double norm = dot_func(grad_saved.data(), grad_saved.data(), dim);
        norm = sqrt(norm);
        for (int i = 0; i < dim; i++)
        {
            grad_saved[i] /= norm;
        }
    }

    Ions_Move_Basic::check_converged(ucell, grad.data(), ofs);
    if (Ions_Move_Basic::converged)
    {
        Ions_Move_Basic::terminate(ucell, ofs);
    }
    else
    {
        this->cal_tradius_sd();
        for (int i = 0; i < dim; i++)
        {
            move[i] = -grad_saved[i] * trust_radius;
        }
        move_atoms(ucell, move.data(), pos_saved.data(), ofs);
        Ions_Move_Basic::update_iter++;
    }

    return;
}

void Ions_Move_SD::cal_tradius_sd() const
{
    static int accepted_number = 0;

    if (Ions_Move_Basic::istep == 1)
    {
        Ions_Move_Basic::trust_radius = Ions_Move_Basic::relax_bfgs_init;
    }
    else if (Ions_Move_Basic::istep > 1)
    {
        if (Ions_Move_Basic::ediff < 0.0)
        {
            accepted_number++;
            if (accepted_number > 3 && accepted_number % 3 == 1)
            {
                Ions_Move_Basic::trust_radius *= 1.5;
            }
        }
        else if (Ions_Move_Basic::ediff >= 0.0) // == 0 means no accept!
        {
            accepted_number = 0;
            Ions_Move_Basic::trust_radius *= 0.5;
        }
    }
    else
    {
        ModuleBase::WARNING_QUIT("Ions_Move_SD::cal_tradius_sd", "istep < 1!");
    }
    if (PARAM.inp.out_level == "ie")
    {
        std::cout << " SD RADIUS (Bohr)     : " << trust_radius << std::endl;
    }
    return;
}
