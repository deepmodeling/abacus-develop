#ifndef GRID_H
#define GRID_H

#include "module_cell/unitcell.h"
#include "sltk_atom.h"
#include "sltk_util.h"

#include <functional>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

typedef std::vector<FAtom> AtomMap;

class Grid
{
  public:
    // Constructors and destructor
    // Grid is Global class,so init it with constant number
    Grid() : test_grid(0){};
    Grid(const int& test_grid_in);
    virtual ~Grid();

    Grid& operator=(Grid&&) = default;

    void init(std::ofstream& ofs, const UnitCell& ucell, const double radius_in, const bool boundary = true);

    // Data
    bool pbc=false; // When pbc is set to false, periodic boundary conditions are explicitly ignored.
    double sradius2=0.0; // searching radius squared (unit:lat0)
    double sradius=0.0;  // searching radius (unit:lat0)
    
    // coordinate range of the input atom (unit:lat0)
    double x_min=0.0;
    double y_min=0.0;
    double z_min=0.0;
    double x_max=0.0;
    double y_max=0.0;
    double z_max=0.0;

    // The algorithm for searching neighboring atoms uses a "box" partitioning method. 
    // Each box has an edge length of sradius, and the number of boxes in each direction is recorded here.
    double box_edge_length=0.0;
    int box_nx=0;
    int box_ny=0;
    int box_nz=0;

    void getBox(int& bx, int& by, int& bz, const double& x, const double& y, const double& z)
    {
        bx = std::floor((x - x_min) / box_edge_length);
        by = std::floor((y - y_min) / box_edge_length);
        bz = std::floor((z - z_min) / box_edge_length);
    }
    // Stores the atoms after box partitioning.
    std::vector<std::vector<std::vector<AtomMap>>> atoms_in_box;

    // Stores the adjacent information of atoms. [ntype][natom][adj list]
    std::vector<std::vector< std::vector<FAtom *> >> all_adj_info;
    void clear_atoms()
    {
        // we have to clear the all_adj_info
        // because the pointers point to the memory in vector atoms_in_box
        all_adj_info.clear();

        atoms_in_box.clear();
    }
    void clear_adj_info()
    {
        // here dont need to free the memory, 
        // because the pointers point to the memory in vector atoms_in_box
        all_adj_info.clear();
    }
    int getGlayerX() const
    {
        return glayerX;
    }
    int getGlayerY() const
    {
        return glayerY;
    }
    int getGlayerZ() const
    {
        return glayerZ;
    }
    int getGlayerX_minus() const
    {
        return glayerX_minus;
    }
    int getGlayerY_minus() const
    {
        return glayerY_minus;
    }
    int getGlayerZ_minus() const
    {
        return glayerZ_minus;
    }
  private:
    int test_grid;

    void setMemberVariables(std::ofstream& ofs_in, const UnitCell& ucell);

    void Construct_Adjacent(const UnitCell& ucell);
    void Construct_Adjacent_near_box(const FAtom& fatom);
    void Construct_Adjacent_final(const FAtom& fatom1, FAtom* fatom2);

    void Check_Expand_Condition(const UnitCell& ucell);
    int glayerX=0;
    int glayerX_minus=0;
    int glayerY=0;
    int glayerY_minus=0;
    int glayerZ=0;
    int glayerZ_minus=0;
};

#endif
