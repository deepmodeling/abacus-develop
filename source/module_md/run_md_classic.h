#ifndef RUN_MD_CLASSIC_H
#define RUN_MD_CLASSIC_H

#include "../module_cell/unitcell_pseudo.h"
#include "../module_esolver/esolver.h"

class Run_MD_CLASSIC
{
public:
    Run_MD_CLASSIC();
    ~Run_MD_CLASSIC();

    void classic_md_line(UnitCell_pseudo &unit_in, ModuleESolver::ESolver *p_esolver);

};

#endif