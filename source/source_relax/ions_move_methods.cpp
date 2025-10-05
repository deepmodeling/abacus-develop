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

void Ions_Move_Methods::allocate(const int &natom)
{
    Ions_Move_Basic::dim = natom * 3;

    if (Ions_Move_Basic::relax_method[0] == "bfgs"&&Ions_Move_Basic::relax_method[1] != "1")
    {
        this->bfgs.allocate();
    }
    else if (Ions_Move_Basic::relax_method[0] == "sd")
    {
        this->sd.allocate();
    }
    else if (Ions_Move_Basic::relax_method[0] == "cg")
    {
        this->cg.allocate();
    }
    else if (Ions_Move_Basic::relax_method[0] == "cg_bfgs")
    {
        this->cg.allocate();
        this->bfgs.allocate(); // added by pengfei  13-8-8
    }
    else if(Ions_Move_Basic::relax_method[0] == "bfgs"&&Ions_Move_Basic::relax_method[1] == "1")
    {
        this->bfgs_trad.allocate(natom);       
    }
    else if(Ions_Move_Basic::relax_method[0] == "lbfgs")
    {
        this->lbfgs.allocate(natom);       
    }
    else
    {
        ModuleBase::WARNING("Ions_Move_Methods::init", "the parameter Ions_Move_Basic::relax_method is not correct.");
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
    if (Ions_Move_Basic::relax_method[0] == "bfgs"&&Ions_Move_Basic::relax_method[1] != "1")
    {
        // move_ions
        // output tau
        // check all symmery
        bfgs.start(ucell, f, etot);
    }
    else if (Ions_Move_Basic::relax_method[0] == "sd")
    {
        sd.start(ucell, f, etot);
    }
    else if (Ions_Move_Basic::relax_method[0] == "cg")
    {
        cg.start(ucell, f, etot);
    }
    else if (Ions_Move_Basic::relax_method[0] == "cg_bfgs")
    {
        cg.start(ucell, f, etot); // added by pengfei 13-8-10
    }
    else if(Ions_Move_Basic::relax_method[0] == "bfgs"&&Ions_Move_Basic::relax_method[1] == "1")
    {
        bfgs_trad.relax_step(f,ucell);        
    }
    else if(Ions_Move_Basic::relax_method[0] == "lbfgs")
    {
        lbfgs.relax_step(f,ucell,etot);        
    }
    else
    {
        ModuleBase::WARNING("Ions_Move_Methods::init", "the parameter Ions_Move_Basic::relax_method is not correct.");
    }
    return;
}
