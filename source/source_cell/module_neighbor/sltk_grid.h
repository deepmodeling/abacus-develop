#ifndef GRID_H
#define GRID_H

#include "source_cell/unitcell.h"
#include "sltk_atom.h"

#include <cmath>
#include <functional>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <vector>

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

    /**
     * @brief Build the expanded-cell neighbor table for one UnitCell.
     *
     * @param [out] ofs Output stream used by ABACUS running logs.
     * @param [in] ucell Unit cell whose atom positions are in lat0-scaled
     * Cartesian coordinates.
     * @param [in] radius_in Search radius in lat0 units.
     * @param [in] boundary Whether periodic boundary images are included.
     *
     * @details The current implementation first estimates periodic image
     * layers, then builds atom images in a Cartesian cell list, and finally
     * constructs all_adj_info for every center atom.
     */
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

    // Cartesian cell-list partitioning. Each bin edge is at least sradius so a
    // 27-bin stencil is sufficient for Euclidean cutoff search.
    double box_edge_length=0.0;
    int box_nx=0;
    int box_ny=0;
    int box_nz=0;

    void getBox(int& bx, int& by, int& bz, const double& x, const double& y, const double& z) const
    {
        if (box_edge_length <= 0.0 || box_nx <= 0 || box_ny <= 0 || box_nz <= 0)
        {
            bx = 0;
            by = 0;
            bz = 0;
            return;
        }

        bx = static_cast<int>(std::floor((x - x_min) / box_edge_length));
        by = static_cast<int>(std::floor((y - y_min) / box_edge_length));
        bz = static_cast<int>(std::floor((z - z_min) / box_edge_length));

        if (bx < 0)
        {
            bx = 0;
        }
        else if (bx >= box_nx)
        {
            bx = box_nx - 1;
        }

        if (by < 0)
        {
            by = 0;
        }
        else if (by >= box_ny)
        {
            by = box_ny - 1;
        }

        if (bz < 0)
        {
            bz = 0;
        }
        else if (bz >= box_nz)
        {
            bz = box_nz - 1;
        }
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

    /**
     * @brief Build expanded atom images and resize neighbor storage.
     */
    void setMemberVariables(std::ofstream& ofs_in, const UnitCell& ucell);

    /**
     * @brief Construct neighbor pointer lists for all atoms in the original cell.
     */
    void Construct_Adjacent(const UnitCell& ucell);

    /**
     * @brief Scan the Cartesian cell-list stencil for one center atom.
     */
    void Construct_Adjacent_near_box(const FAtom& fatom, std::vector<FAtom*>& adjacent_atoms);

    /**
     * @brief Append fatom2 to adjacent_atoms when the current distance rule
     * accepts the pair.
     */
    void Construct_Adjacent_final(const FAtom& fatom1,
                                  FAtom* fatom2,
                                  std::vector<FAtom*>& adjacent_atoms) const;

    /**
     * @brief Estimate positive and negative periodic image layers needed for sradius.
     */
    void Check_Expand_Condition(const UnitCell& ucell);
    int glayerX=0;
    int glayerX_minus=0;
    int glayerY=0;
    int glayerY_minus=0;
    int glayerZ=0;
    int glayerZ_minus=0;
};

#endif
