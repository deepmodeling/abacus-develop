#ifndef GRID_MESHBALL_H
#define GRID_MESHBALL_H

#include "grid_bigcell.h"

class Grid_MeshBall : public Grid_BigCell
{
  public:
    // the radius of meshball.
    double meshball_radius;

    // number of meshcells in meshball.
    // generally, meshball_radius = orbital_rmax;
    int meshball_ncells;

    // cartesian coordinates of meshball.
    double** meshball_positions;

    bool flag_mp;

  protected:
    Grid_MeshBall();
    ~Grid_MeshBall();

    // used in index2normal
    int* index_ball;

    // init the meshball radius,
    // search each meshcell of this meshball.
    void init_meshball(void);

    void delete_meshball_positions(void); // LiuXh add 2018-12-14

  private:
    double deal_with_atom_spillage(const double* pos);
};
#endif
