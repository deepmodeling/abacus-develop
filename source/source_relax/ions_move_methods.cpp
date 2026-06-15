#include "ions_move_methods.h"

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
