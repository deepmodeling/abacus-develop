#include "ions_move_methods.h"

#include "ions_move_basic.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"

Ions_Move_Methods::Ions_Move_Methods()
{
}
Ions_Move_Methods::~Ions_Move_Methods()
{
}

void Ions_Move_Methods::allocate(const int &natom)
{
    Ions_Move_Basic::dim = natom * 3;

    if (GlobalV::RELAX_METHOD == "bfgs")
    {
        this->bfgs.allocate();
    }
    else if (GlobalV::RELAX_METHOD == "sd")
    {
        this->sd.allocate();
    }
    else if (GlobalV::RELAX_METHOD == "cg")
    {
        this->cg.allocate();
    }
    else if (GlobalV::RELAX_METHOD == "cg_bfgs")
    {
        this->cg.allocate();
        this->bfgs.allocate(); // added by pengfei  13-8-8
    }
    else
    {
        ModuleBase::WARNING("Ions_Move_Methods::init", "the parameter GlobalV::RELAX_METHOD is not correct.");
    }
    return;
}

// void Ions_Move_Methods::cal_movement(const int &istep, const ModuleBase::matrix &f, const double &etot)
void Ions_Move_Methods::cal_movement(const int &istep,
                                     const int &force_step,
                                     const ModuleBase::matrix &f,
                                     const double &etot,
                                     UnitCell &ucell)
{
    ModuleBase::TITLE("Ions_Move_Methods", "init");

    // Ions_Move_Basic::istep = istep;
    Ions_Move_Basic::istep = force_step;

    if (GlobalV::RELAX_METHOD == "bfgs")
    {
        // move_ions
        // output tau
        // check all symmery
        bfgs.start(ucell, f, etot);
    }
    else if (GlobalV::RELAX_METHOD == "sd")
    {
        sd.start(ucell, f, etot);
    }
    else if (GlobalV::RELAX_METHOD == "cg")
    {
        cg.start(ucell, f, etot);
    }
    else if (GlobalV::RELAX_METHOD == "cg_bfgs")
    {
        cg.start(ucell, f, etot); // added by pengfei 13-8-10
    }
    else
    {
        ModuleBase::WARNING("Ions_Move_Methods::init", "the parameter GlobalV::RELAX_METHOD is not correct.");
    }
    return;
}
