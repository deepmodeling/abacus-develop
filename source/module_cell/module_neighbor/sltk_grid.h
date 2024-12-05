#ifndef GRID_H
#define GRID_H

#include "module_cell/unitcell.h"
#include "sltk_atom.h"
#include "sltk_util.h"

#include <functional>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

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

    std::vector<std::vector<std::vector<std::vector<FAtom>>>> Cell; // dx , dy ,dz is cell number in each direction,respectly.
    std::vector<std::vector<std::vector<FAtom *>>> all_adj_info;

    int getGlayerX()
    {
        return glayerX;
    }
    int getGlayerY()
    {
        return glayerY;
    }
    int getGlayerZ()
    {
        return glayerZ;
    }
    int getGlayerX_minus()
    {
        return glayerX_minus;
    }
    int getGlayerY_minus()
    {
        return glayerY_minus;
    }
    int getGlayerZ_minus()
    {
        return glayerZ_minus;
    }
    

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

  private:
    const int test_grid;

    void setMemberVariables(std::ofstream& ofs_in);
    void Build_Hash_Table(const UnitCell& ucell);

    void Construct_Adjacent_expand(const UnitCell & ucell);

    void Construct_Adjacent_expand_periodic(const FAtom& fatom);

    void Construct_Adjacent_final(const FAtom& fatom1,
                                  FAtom& fatom2);
};

#endif
