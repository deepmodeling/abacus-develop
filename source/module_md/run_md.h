#ifndef RUN_MD_CLASSIC_H
#define RUN_MD_CLASSIC_H

#include "md_para.h"
#include "module_esolver/esolver.h"

class Run_MD
{
  public:
    Run_MD();
    ~Run_MD();

    void md_line(UnitCell &unit_in, ModuleESolver::ESolver *p_esolver, MD_parameters &md_para);
};

#endif