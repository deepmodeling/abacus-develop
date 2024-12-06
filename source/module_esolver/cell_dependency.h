#ifndef CELL_DEPENDENCY
#define CELL_DEPENDENCY

#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_hamilt_lcao/hamilt_lcaodft/record_adj.h"
#include "module_lr/utils/lr_util.hpp"

#include <memory>

class Cell_Dependency
{
  public:
    Cell_Dependency() {};
    ~Cell_Dependency() {};

    Record_adj ra;
    Grid_Driver grid_d;
};

#endif // CELL_DEPENDENCY
