//==========================================================
// AUTHOR : mohan
// DATE : 2021-02-01
//==========================================================
#ifndef VARIABLE_CELL_H
#define VARIABLE_CELL_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../input.h"
#include "module_esolver/esolver.h"

class Variable_Cell
{

	public:

    Variable_Cell();
    ~Variable_Cell();

    static void init_after_vc(); //LiuXh add 20180515

};

#endif
