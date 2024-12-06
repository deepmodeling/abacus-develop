#ifndef ATOM_INPUT_H
#define ATOM_INPUT_H

#include "module_cell/unitcell.h"
#include "sltk_atom.h"

class Atom_input
{
  public:
    //==========================================================
    // Constructors and destructor
    //==========================================================
    Atom_input(std::ofstream& ofs_in,
               const UnitCell& ucell,
               const int amount = 0,       // number of atoms
               const int ntype = 0,        // number of atom_types
               const bool boundary = true, // 1 : periodic ocndition
               const double radius_in = 0, // searching radius
               const int& test_atom_in = 0 // caoyu reconst 2021-05-24
    );
    ~Atom_input();

  public:
    bool getExpandFlag() const
    {
        return expand_flag;
    }

    int getBoundary() const
    {
        return periodic_boundary;
    }

    double getRadius() const
    {
        return radius;
    }

    //==========================================================
    //
    //==========================================================
    double minX() const
    {
		return x_min;
    }

    double minY() const
    {
		return y_min;
    }

    double minZ() const
    {
		return z_min;
    }

    //==========================================================
    //
    //==========================================================
    int getCell_nX() const
    {
        return cell_nx;
    }

    int getCell_nY() const
    {
        return cell_ny;
    }

    int getCell_nZ() const
    {
        return cell_nz;
    }

    //==========================================================
    //
    //==========================================================
    int getGrid_layerX() const
    {
        return glayerX;
    }

    int getGrid_layerX_minus() const
    {
        return glayerX_minus;
    }

    int getGrid_layerY() const
    {
        return glayerY;
    }

    int getGrid_layerY_minus() const
    {
        return glayerY_minus;
    }

    int getGrid_layerZ() const
    {
        return glayerZ;
    }

    int getGrid_layerZ_minus() const
    {
        return glayerZ_minus;
    }

  private:
    int test_atom_input; // caoyu reconst 2021-05-24
    bool periodic_boundary;

    double radius;

    double x_min;
    double y_min;
    double z_min;
    double x_max;
    double y_max;
    double z_max;
    //==========================================================
    // MEMBRE FUNCTION :
    // NAME : Check_Expand_Condition
    //==========================================================
    void Check_Expand_Condition(const UnitCell& ucell);
    bool expand_flag;
    int glayerX;
    int glayerX_minus;
    int glayerY;
    int glayerY_minus;
    int glayerZ;
    int glayerZ_minus;

    //==========================================================
    // MEMBRE FUNCTION :
    // NAME : Expand_Grid
    //==========================================================
    void calculate_cells();
    int cell_nx;
    int cell_ny;
    int cell_nz;
};

#endif
