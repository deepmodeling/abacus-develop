#ifndef LCAO_gen_fixedH_H
#define LCAO_gen_fixedH_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_basis/module_ao/ORB_gen_tables.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "module_basis/module_ao/ORB_gen_tables.h"
#include "module_base/vector3.h"

class LCAO_gen_fixedH
{

  public:
    LCAO_Matrix* LM;

    LCAO_gen_fixedH();
    ~LCAO_gen_fixedH();

};

#endif
