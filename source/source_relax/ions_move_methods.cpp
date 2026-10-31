#include "ions_move_methods.h"

#include "ions_move_basic.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_io/module_parameter/parameter.h"

#include <algorithm>

Ions_Move_Methods::Ions_Move_Methods()
{
}
Ions_Move_Methods::~Ions_Move_Methods()
{
}

void Ions_Move_Methods::allocate(const int& natom, const std::string& relax_method)
{
    if (natom <= 0)
    {
        ModuleBase::WARNING_QUIT("Ions_Move_Methods::allocate", "natom must be greater than 0.");
    }
    Ions_Move_Basic::dim = natom * 3;

    if (relax_method == "bfgs")
    {
        this->bfgs.allocate();
    }
    else if (relax_method == "sd")
    {
        this->sd.allocate();
    }
    else if (relax_method == "cg_trad")
    {
        this->cg.allocate(Ions_Move_Basic::dim);
    }
    else if (relax_method == "cg_bfgs")
    {
        this->cg.allocate(Ions_Move_Basic::dim);
        this->bfgs.allocate();
    }
    else if (relax_method == "bfgs_trad")
    {
        this->bfgs_trad.allocate(natom);
    }
    else if (relax_method == "lbfgs")
    {
        this->lbfgs.allocate(natom, PARAM.inp.relax_bfgs_rmax, PARAM.inp.out_level);
    }
    else
    {
        ModuleBase::WARNING("Ions_Move_Methods::init", "the parameter relax_method is not correct.");
    }
    return;
}

// void Ions_Move_Methods::cal_movement(const int &istep, const ModuleBase::matrix &f, const double &etot)
void Ions_Move_Methods::cal_movement(const int& istep,
                                     const int& force_step,
                                     const ModuleBase::matrix& f,
                                     const double& etot,
                                     UnitCell& ucell,
                                     std::ofstream& ofs,
                                     std::string& relax_method,
                                     const Relax_Criteria& criteria)
{
    ModuleBase::TITLE("Ions_Move_Methods", "init");
    if (relax_method == "bfgs")
    {
        converged_ = bfgs.start(ucell, f, etot, force_step, update_iter_, ofs, etot_info_, criteria);
    }
    else if (relax_method == "sd")
    {
        converged_ = sd.start(ucell, f, etot, force_step, update_iter_, ofs, etot_info_, criteria);
    }
    else if (relax_method == "cg_trad")
    {
        converged_ = cg.start(ucell, f, etot, force_step, update_iter_, ofs, etot_info_, relax_method, criteria);
    }
    else if (relax_method == "cg_bfgs")
    {
        converged_ = cg.start(ucell, f, etot, force_step, update_iter_, ofs, etot_info_, relax_method, criteria);
    }
    else if (relax_method == "bfgs_trad")
    {
        converged_ = bfgs_trad.relax_step(f, ucell, ofs, criteria);
    }
    else if (relax_method == "lbfgs")
    {
        converged_ = lbfgs.relax_step(f, ucell, etot, ofs, criteria);
    }
    else
    {
        ModuleBase::WARNING("Ions_Move_Methods::init", "the parameter relax_method is not correct.");
        converged_ = false;
    }
    return;
}

void Ions_Move_Methods::reset_after_cell_change(const std::string& relax_method, std::ofstream& ofs)
{
    ModuleBase::TITLE("Ions_Move_Methods", "reset_after_cell_change");

    if (relax_method.empty())
    {
        return;
    }

    const auto reset_common_state = [this]() {
        this->converged_ = false;
        this->update_iter_ = 0;
        std::fill(this->etot_info_.begin(), this->etot_info_.end(), 0.0);
    };

    if (relax_method == "bfgs")
    {
        reset_common_state();
        this->bfgs.reset();
        ofs << " Reset ionic BFGS history after cell change." << std::endl;
    }
    else if (relax_method == "bfgs_trad")
    {
        reset_common_state();
        this->bfgs_trad.reset();
        ofs << " Reset traditional ionic BFGS history after cell change." << std::endl;
    }
}
