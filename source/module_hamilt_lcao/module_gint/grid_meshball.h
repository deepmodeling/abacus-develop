#ifndef GRID_MESHBALL_H
#define GRID_MESHBALL_H

#include "grid_bigcell.h"

class Grid_MeshBall : public Grid_BigCell
{
  public:
    Grid_MeshBall();
    ~Grid_MeshBall();
    // cartesian coordinates of meshball.
    std::vector<std::vector<double>> meshball_positions;

  protected:
    // number of meshcells in meshball.
    int meshball_ncells;
    // used in index2normal
    std::vector<int> index_ball;
    // search each meshcell of this meshball.
    void init_meshball(void);

  private:
    // flag_mp determines whether meshball_positions is configured.
    bool flag_mp;
    // init the meshball radius.
    double meshball_radius;
    // Handle as a truncation function.
    double deal_with_atom_spillage(const double* pos);
};
#endif
