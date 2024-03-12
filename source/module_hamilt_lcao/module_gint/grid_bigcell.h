#ifndef GRID_BIGCELL_H
#define GRID_BIGCELL_H

#include "grid_meshcell.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/matrix3.h"

class Grid_BigCell : public Grid_MeshCell {
public:
  Grid_BigCell();
  ~Grid_BigCell();

  // save the relative cartesian position
  // to bigcell of each atom.
  double **tau_in_bigcell;

protected:
  //---------------------------------
  void init_big_latvec(void);

  double bigcell_vec1[3];
  double bigcell_vec2[3];
  double bigcell_vec3[3];

  ModuleBase::Matrix3 bigcell_latvec0;
  ModuleBase::Matrix3 bigcell_GT;
  //---------------------------------

  //---------------------------------
  void init_grid_expansion(void);

  // get the max radius of all orbitals.
  // which will use to generate grid expansion,
  // and  the meshball.
  double orbital_rmax;

  // Measure the distance between two bigcell
  // in each direction.
  double bigcell_dx;
  double bigcell_dy;
  double bigcell_dz;

  // the added number of bigcelli each direction.
  int dxe;
  int dye;
  int dze;

  // expansion grid dimension.
  int nxe;
  int nye;
  int nze;
  int nxyze;
  //---------------------------------

  //---------------------------------
  void init_tau_in_bigcell(void);

  // this flag will be false at first and turned to true after memory of
  // tau_in_meshcell has been allocated.
  bool flag_tib;

  int *index_atom;
  //---------------------------------

  //---------------------------------
  void grid_expansion_index(bool f2normal, int *target) const;
  //---------------------------------
};
#endif
