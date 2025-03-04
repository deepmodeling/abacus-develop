#ifndef ELECSATE_PRINT_H
#define ELECSATE_PRINT_H
#include "module_elecstate/elecstate.h"
namespace elecstate
{
    void  print_band(const ModuleBase::matrix& ekb, 
                const ModuleBase::matrix& wg,
                const K_Vectors* klist,
                const int& ik, 
                const int& printe, 
                const int& iter);
                
}
#endif