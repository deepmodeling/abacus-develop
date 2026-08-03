/**
 * @file sltk_grid.h
 * @brief Grid class for neighbor search.
 */
#ifndef GRID_H
#define GRID_H

#include "source_cell/unitcell.h"
#include "sltk_atom.h"

#include <functional>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

typedef std::vector<FAtom> AtomMap;

/**
 * @brief Grid class for neighbor search.
 *
 * The algorithm for searching neighboring atoms uses a "box" partitioning method.
 * Each box has an edge length of sradius, and the number of boxes in each direction is recorded.
 */
class Grid
{
  public:
    /**
     * @brief Default constructor.
     *
     * Grid is Global class, so init it with constant number.
     */
    Grid() : test_grid(0){};

    /**
     * @brief Constructor with test flag.
     *
     * @param test_grid_in test flag
     */
    Grid(const int& test_grid_in);

    /**
     * @brief Destructor.
     */
    virtual ~Grid();

    Grid& operator=(Grid&&) = default;

    /**
     * @brief Initialize the grid.
     *
     * @param ofs output file stream
     * @param ucell unit cell
     * @param radius_in searching radius
     * @param boundary whether to apply boundary conditions
     */
    void init(std::ofstream& ofs, const UnitCell& ucell, const double radius_in, const bool boundary = true);

    /// @brief Data

    /// @brief When pbc is set to false, periodic boundary conditions are explicitly ignored.
    bool pbc=false;

    /// @brief searching radius squared (unit:lat0)
    double sradius2=0.0;

    /// @brief searching radius (unit:lat0)
    double sradius=0.0;
    
    /// @brief coordinate range of the input atom (unit:lat0)
    double x_min=0.0;
    double y_min=0.0;
    double z_min=0.0;
    double x_max=0.0;
    double y_max=0.0;
    double z_max=0.0;

    /// @brief box edge length (equal to sradius)
    double box_edge_length=0.0;

    /// @brief number of boxes in x direction
    int box_nx=0;

    /// @brief number of boxes in y direction
    int box_ny=0;

    /// @brief number of boxes in z direction
    int box_nz=0;

    /**
     * @brief Get box indices for given coordinates.
     *
     * @param bx box index in x direction (output)
     * @param by box index in y direction (output)
     * @param bz box index in z direction (output)
     * @param x x coordinate
     * @param y y coordinate
     * @param z z coordinate
     */
    void getBox(int& bx, int& by, int& bz, const double& x, const double& y, const double& z)
    {
        bx = std::floor((x - x_min) / box_edge_length);
        by = std::floor((y - y_min) / box_edge_length);
        bz = std::floor((z - z_min) / box_edge_length);
    }

    /// @brief Stores the atoms after box partitioning.
    std::vector<std::vector<std::vector<AtomMap>>> atoms_in_box;

    /// @brief Stores the adjacent information of atoms. [ntype][natom][adj list]
    std::vector<std::vector< std::vector<FAtom *> >> all_adj_info;

    /**
     * @brief Clear all atoms and adjacent information.
     *
     * We have to clear the all_adj_info because the pointers point to the memory in vector atoms_in_box.
     */
    void clear_atoms()
    {
        all_adj_info.clear();
        atoms_in_box.clear();
    }

    /**
     * @brief Clear adjacent information only.
     *
     * Here we don't need to free the memory because the pointers point to the memory in vector atoms_in_box.
     */
    void clear_adj_info()
    {
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
    int test_grid; ///< test flag

    /**
     * @brief Set member variables.
     *
     * @param ofs_in output file stream
     * @param ucell unit cell
     */
    void setMemberVariables(std::ofstream& ofs_in, const UnitCell& ucell);

    /**
     * @brief Construct adjacent atom information.
     *
     * @param ucell unit cell
     */
    void Construct_Adjacent(const UnitCell& ucell);

    /**
     * @brief Construct adjacent atom information for nearby boxes.
     *
     * @param fatom atom for which to find neighbors
     */
    void Construct_Adjacent_near_box(const FAtom& fatom);

    /**
     * @brief Finalize adjacent atom information for a pair of atoms.
     *
     * @param fatom1 first atom
     * @param fatom2 second atom
     */
    void Construct_Adjacent_final(const FAtom& fatom1, FAtom* fatom2);

    /**
     * @brief Check expansion condition for periodic images.
     *
     * @param ucell unit cell
     */
    void Check_Expand_Condition(const UnitCell& ucell);

    int glayerX=0;       ///< number of periodic images in positive x direction
    int glayerX_minus=0; ///< number of periodic images in negative x direction
    int glayerY=0;       ///< number of periodic images in positive y direction
    int glayerY_minus=0; ///< number of periodic images in negative y direction
    int glayerZ=0;       ///< number of periodic images in positive z direction
    int glayerZ_minus=0; ///< number of periodic images in negative z direction
};

#endif
