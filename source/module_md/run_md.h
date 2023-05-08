#ifndef RUN_MD_H
#define RUN_MD_H

#include "md_para.h"
#include "module_esolver/esolver.h"

namespace Run_MD
{
void md_line(UnitCell& unit_in, ModuleESolver::ESolver* p_esolver, MD_parameters& md_para);
} // namespace Run_MD

#endif