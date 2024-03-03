#ifndef GRID_MESHK_H
#define GRID_MESHK_H
#include "module_base/global_function.h"
#include "module_base/global_variable.h"

class Grid_MeshK {
public:
  // the max and the min unitcell.
  int maxu1, maxu2, maxu3;
  int minu1, minu2, minu3;

  // the number of unitcells.
  int nu1, nu2, nu3;
  int nutot;

  // from 1D index to unitcell.
  int *ucell_index2x;
  int *ucell_index2y;
  int *ucell_index2z;

  int cal_Rindex(const int &u1, const int &u2, const int &u3) const;

protected:
  Grid_MeshK();
  ~Grid_MeshK();
};

#endif
