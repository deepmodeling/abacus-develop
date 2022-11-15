#ifndef RUN_MD_CLASSIC_H
#define RUN_MD_CLASSIC_H

#include "../module_cell/unitcell.h"
#include "../module_esolver/esolver.h"

class Run_MD
{
public:
    Run_MD();
    ~Run_MD();

    void md_line(UnitCell &unit_in, ModuleESolver::ESolver *p_esolver);

};

#endif