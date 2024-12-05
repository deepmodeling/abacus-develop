#ifndef GRID_H
#define GRID_H

#include "module_cell/unitcell.h"
#include "sltk_atom.h"
#include "sltk_util.h"

#include <functional>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

typedef std::vector<std::vector<FAtom>> AtomMap;

struct CellSet
{
    AtomMap atom_map;
    int in_grid[3];
    CellSet();
};

//==========================================================
// CLASS NAME :
// Atom_input : defined elsewhere
//==========================================================

//==========================================================
// CLASS NAME :
// Grid :
//==========================================================

class Grid
{
  public:
    // Constructors and destructor
    // Grid is Global class,so init it with constant number
    Grid() : test_grid(0){};
    Grid(const int& test_grid_in);
    virtual ~Grid();

void init(std::ofstream& ofs_in, const UnitCell& ucell, const double radius, bool pbc_in);

    // Data
    bool pbc; // periodic boundary condition

    double sradius2; // searching radius squared
    double sradius;  // searching radius

    int cell_nx;
    int cell_ny;
    int cell_nz;

    int true_cell_x;
    int true_cell_y;
    int true_cell_z;
    double x_min;
    double y_min;
    double z_min;
    double x_max;
    double y_max;
    double z_max;

    void Check_Expand_Condition(const UnitCell& ucell);
    int glayerX;
    int glayerX_minus;
    int glayerY;
    int glayerY_minus;
    int glayerZ;
    int glayerZ_minus;

    std::vector<std::vector<std::vector<CellSet>>> Cell; // dx , dy ,dz is cell number in each direction,respectly.

    int getCellX() const
    {
        return cell_nx;
    }
    int getCellY() const
    {
        return cell_ny;
    }
    int getCellZ() const
    {
        return cell_nz;
    }
    int getTrueCellX() const
    {
        return true_cell_x;
    }
    int getTrueCellY() const
    {
        return true_cell_y;
    }
    int getTrueCellZ() const
    {
        return true_cell_z;
    }

  private:
    const int test_grid;

    void setMemberVariables(std::ofstream& ofs_in);

    void setBoundaryAdjacent(std::ofstream& ofs_in);

    //==========================================================
    void Build_Hash_Table(const UnitCell& ucell);

    //==========================================================

    void Construct_Adjacent_expand(const int i, const int j, const int k);

    void Construct_Adjacent_expand_periodic(const int i, const int j, const int k, FAtom& fatom);

    void Construct_Adjacent_begin();
    void Construct_Adjacent_nature(const int i, const int j, const int k, FAtom& fatom1);
    void Construct_Adjacent_periodic(const int i, const int j, const int k, FAtom& fatom1);
    void Construct_Adjacent_final(const int i,
                                  const int j,
                                  const int k,
                                  FAtom& fatom1,
                                  const int i2,
                                  const int j2,
                                  const int k2,
                                  FAtom& fatom2);
};

#endif
