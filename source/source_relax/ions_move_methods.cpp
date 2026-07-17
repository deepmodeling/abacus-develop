#include "ions_move_methods.h"

#include <algorithm>

#include "ions_move_basic.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"


Ions_Move_Methods::Ions_Move_Methods()
{
}
Ions_Move_Methods::~Ions_Move_Methods()
{
}

void Ions_Move_Methods::allocate(const int &natom, const std::string& relax_method_0, const std::string& relax_method_1)
{
    if (natom <= 0)
    {
        ModuleBase::WARNING_QUIT("Ions_Move_Methods::allocate", "natom must be greater than 0.");
    }
    Ions_Move_Basic::dim = natom * 3;

    if (relax_method_0 == "bfgs" && relax_method_1 != "1")
    {
        this->bfgs.allocate();
    }
    else if (relax_method_0 == "sd")
    {
        this->sd.allocate();
    }
    else if (relax_method_0 == "cg")
    {
        this->cg.allocate(Ions_Move_Basic::dim);
    }
    else if (relax_method_0 == "cg_bfgs")
    {
        this->cg.allocate(Ions_Move_Basic::dim);
        this->bfgs.allocate();
    }
    else if(relax_method_0 == "bfgs" && relax_method_1 == "1")
    {
        this->bfgs_trad.allocate(natom);       
    }
    else if(relax_method_0 == "lbfgs")
    {
        this->lbfgs.allocate(natom);       
    }
    else
    {
        ModuleBase::WARNING("Ions_Move_Methods::init", "the parameter relax_method is not correct.");
    }
    return;
}

// void Ions_Move_Methods::cal_movement(const int &istep, const ModuleBase::matrix &f, const double &etot)
void Ions_Move_Methods::cal_movement(const int &istep,
                                     const int &force_step,
                                     const ModuleBase::matrix &f,
                                     const double &etot,
                                     UnitCell &ucell,
                                     std::ofstream& ofs,
                                     std::vector<std::string>& relax_method)
{
    ModuleBase::TITLE("Ions_Move_Methods", "init");
    if (relax_method[0] == "bfgs" && relax_method[1] != "1")
    {
        converged_ = bfgs.start(ucell, f, etot, force_step, update_iter_, ofs, etot_info_);
    }
    else if (relax_method[0] == "sd")
    {
        converged_ = sd.start(ucell, f, etot, force_step, update_iter_, ofs, etot_info_);
    }
    else if (relax_method[0] == "cg")
    {
        converged_ = cg.start(ucell, f, etot, force_step, update_iter_, ofs, etot_info_, relax_method);
    }
    else if (relax_method[0] == "cg_bfgs")
    {
        converged_ = cg.start(ucell, f, etot, force_step, update_iter_, ofs, etot_info_, relax_method);
    }
    else if (relax_method[0] == "bfgs" && relax_method[1] == "1")
    {
        converged_ = bfgs_trad.relax_step(f, ucell, ofs);        
    }
    else if (relax_method[0] == "lbfgs")
    {
        converged_ = lbfgs.relax_step(f, ucell, etot, ofs);        
    }
    else
    {
        ModuleBase::WARNING("Ions_Move_Methods::init", "the parameter relax_method is not correct.");
        converged_ = false;
    }
    return;
}

void Ions_Move_Methods::reset_after_cell_change(const std::vector<std::string>& relax_method, std::ofstream& ofs)
{
    ModuleBase::TITLE("Ions_Move_Methods", "reset_after_cell_change");

    if (relax_method.empty())
    {
        return;
    }

    const std::string method = relax_method[0];
    const std::string method_arg = relax_method.size() > 1 ? relax_method[1] : "";
    const auto reset_common_state = [this]() {
        this->converged_ = false;
        this->update_iter_ = 0;
        std::fill(this->etot_info_.begin(), this->etot_info_.end(), 0.0);
    };

    if (method == "bfgs" && method_arg != "1")
    {
        reset_common_state();
        this->bfgs.reset();
        ofs << " Reset ionic BFGS history after cell change." << std::endl;
    }
    else if (method == "bfgs" && method_arg == "1")
    {
        reset_common_state();
        this->bfgs_trad.reset();
        ofs << " Reset traditional ionic BFGS history after cell change." << std::endl;
    }
}
