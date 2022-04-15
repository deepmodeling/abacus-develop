#ifndef RUN_MD_CLASSIC_H
#define RUN_MD_CLASSIC_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_cell/unitcell_pseudo.h"

class Run_MD_CLASSIC
{
public:
    Run_MD_CLASSIC();
    ~Run_MD_CLASSIC();

    UnitCell_pseudo ucell_c;
    

    void classic_md_line(void);
    
};

#endif