#ifndef CAL_UX_H
#define CAL_UX_H

#include "module_cell/unitcell.h"

namespace elecstate {

    void cal_ux(UnitCell& ucell);
    
    bool judge_parallel(double a[3], ModuleBase::Vector3<double> b);

}

#endif