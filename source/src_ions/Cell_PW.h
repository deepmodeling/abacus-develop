//================================
// AUTHOR : zhengdy
// DATE : 2021-06-22
//================================
#ifndef CELL_PW_H
#define CELL_PW_H

#include "module_esolver/esolver.h"

class Cell_PW
{
public:
    Cell_PW(){};
    ~Cell_PW(){};
    
    void opt_cells_pw(ModuleESolver::ESolver *p_esolver);
};

#endif