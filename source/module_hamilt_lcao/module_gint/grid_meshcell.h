#ifndef GRID_MESHCELL_H
#define GRID_MESHCELL_H
#include "grid_meshk.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/matrix3.h"

class Grid_MeshCell : public Grid_MeshK {
public:
  // vectors of meshcell.
  double meshcell_vec1[3];
  double meshcell_vec2[3];
  double meshcell_vec3[3];
  ModuleBase::Matrix3 meshcell_latvec0;
  ModuleBase::Matrix3 meshcell_GT;

  double **meshcell_pos;
  bool allocate_pos;

protected:
  Grid_MeshCell();
  ~Grid_MeshCell();

  void set_grid_dim(const int &ncx_in, const int &ncy_in, const int &ncz_in,
                    const int &bx_in, const int &by_in, const int &bz_in,
                    const int &nbx_in, const int &nby_in, const int &nbz_in,
                    const int &nbxx_in, const int &nbzp_start_in,
                    const int &nbzp_in);

  void init_latvec(void);
  void init_meshcell_pos(void);

  void cal_extended_cell(const int &dxe, const int &dye, const int &dze);

public:
  int ncx, ncy, ncz, ncxyz;
  int bx, by, bz, bxyz;
  int nbx, nby, nbz, nbxyz;
  int nbxx;
  int nbzp_start, nbzp;
};

#endif
