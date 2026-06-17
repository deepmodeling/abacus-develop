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

bool Ions_Move_SD::start(UnitCell& ucell, const ModuleBase::matrix& force, const double& etot_in, const int istep, int& update_iter, std::ofstream& ofs, std::vector<double>& etot_info)
{
    ModuleBase::TITLE("Ions_Move_SD", "start");

    assert(dim > 0);
    assert(grad_saved.size() == static_cast<size_t>(dim));
    assert(pos_saved.size() == static_cast<size_t>(dim));

    std::vector<double> pos(dim, 0.0);
    std::vector<double> grad(dim, 0.0);
    std::vector<double> move(dim, 0.0);

    bool judgement = false;
    Ions_Move_Basic::setup_etot(etot_in, judgement, istep, ofs, etot_info);
    Ions_Move_Basic::setup_gradient(ucell, force, pos.data(), grad.data(), ofs);

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
        double norm = dot_func(grad_saved.data(), grad_saved.data(), dim);
        norm = sqrt(norm);
        for (int i = 0; i < dim; i++)
        {
            grad_saved[i] /= norm;
        }
    }

    bool converged = Ions_Move_Basic::check_converged(ucell, grad.data(), update_iter, ofs, etot_info);
    if (converged)
    {
        Ions_Move_Basic::terminate(converged, update_iter, ucell, istep, ofs);
        return true;
    }
    else
    {
        this->cal_tradius_sd(istep, etot_info);
        for (int i = 0; i < dim; i++)
        {
            move[i] = -grad_saved[i] * trust_radius;
        }
        move_atoms(ucell, move.data(), pos_saved.data(), ofs);
        update_iter++;
        return false;
    }
}

void Ions_Move_SD::cal_tradius_sd(const int istep, std::vector<double>& etot_info) const
{
    static int accepted_number = 0;

    if (istep == 1)
    {
        Ions_Move_Basic::trust_radius = Ions_Move_Basic::relax_bfgs_init;
    }
    else if (istep > 1)
    {
        const double ediff = etot_info[0] - etot_info[1];
        if (ediff < 0.0)
        {
            accepted_number++;
            if (accepted_number > 3 && accepted_number % 3 == 1)
            {
                Ions_Move_Basic::trust_radius *= 1.5;
            }
        }
        else if (ediff >= 0.0)
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
