/**
 * @file sltk_grid_driver.h
 * @brief Grid_Driver class for neighbor search interface.
 */
#ifndef GRID_DRIVER_H
#define GRID_DRIVER_H

#include "source_base/global_function.h"
#include "source_base/vector3.h"
#include "source_cell/unitcell.h"
#include "sltk_grid.h"

#include <memory>
#include <stdexcept>
#include <tuple>

/**
 * @brief Struct of array for packing the Adjacent atom information.
 */
class AdjacentAtomInfo
{
  public:
    /**
     * @brief Default constructor.
     */
    AdjacentAtomInfo() : adj_num(0)
    {
    }

    int adj_num;                         ///< number of adjacent atoms
    std::vector<int> ntype;              ///< types of adjacent atoms
    std::vector<int> natom;              ///< indices of adjacent atoms
    std::vector<ModuleBase::Vector3<double>> adjacent_tau; ///< positions of adjacent atoms
    std::vector<ModuleBase::Vector3<int>> box; ///< box indices of adjacent atoms

    /**
     * @brief Clear all adjacent atom information.
     */
    void clear()
    {
        adj_num = 0;
        ntype.clear();
        natom.clear();
        adjacent_tau.clear();
        box.clear();
    }
};

/**
 * @brief Filter adjacent atoms based on boolean mask.
 *
 * @param is_adj boolean mask indicating which atoms are adjacent
 * @param adjs adjacent atom information to filter
 */
void filter_adjs(const std::vector<bool>& is_adj, AdjacentAtomInfo& adjs);

/**
 * @brief Grid_Driver class for neighbor search interface.
 *
 * This class provides the user interface for finding adjacent atoms.
 */
class Grid_Driver : public Grid
{
  public:
    /**
     * @brief Default constructor.
     */
    Grid_Driver(){ test_deconstructor = false; };

    /**
     * @brief Constructor with test flags.
     *
     * @param test_d_in test deconstructor flag
     * @param test_grid_in test grid flag
     */
    Grid_Driver(const int& test_d_in, const int& test_grid_in);

    /**
     * @brief Destructor.
     */
    ~Grid_Driver();

    Grid_Driver& operator=(Grid_Driver&&) = default;

    /**
     * @brief Find adjacent atoms for a given atom.
     *
     * @note This design makes Grid_Driver compatible with multi-thread usage:
     *       1. Find_atom stores results in Grid_Driver::adj_info by default.
     *       2. And stores results into parameter adjs when adjs is NOT NULL.
     *
     * @param ucell unit cell
     * @param ntype atom type
     * @param nnumber atom index within type
     * @param adjs optional output adjacent atom information
     */
    void Find_atom(const UnitCell& ucell,
                   const int ntype,
                   const int nnumber,
                   AdjacentAtomInfo* adjs = nullptr) const;

    /**
     * @brief Find adjacent atoms for a given cartesian position (deprecated).
     *
     * @deprecated This interface is deprecated, please use Find_atom above.
     * @note cartesian_posi and ucell are deprecated 20241204 zhanghaochong
     *
     * @param ucell unit cell
     * @param cartesian_posi cartesian position
     * @param ntype atom type
     * @param nnumber atom index within type
     * @param adjs optional output adjacent atom information
     */
    void Find_atom(const UnitCell& ucell,
                   const ModuleBase::Vector3<double>& cartesian_posi,
                   const int& ntype,
                   const int& nnumber,
                   AdjacentAtomInfo* adjs = nullptr) const;

    /**
     * @brief Get the number of adjacent atoms.
     * @return number of adjacent atoms
     */
    const int& getAdjacentNum() const
    {
        return adj_info.adj_num;
    }

    /**
     * @brief Get the type of an adjacent atom.
     * @param i index of adjacent atom
     * @return type of adjacent atom
     */
    const int& getType(const int i) const
    {
        return adj_info.ntype[i];
    }

    /**
     * @brief Get the index of an adjacent atom.
     * @param i index of adjacent atom
     * @return index of adjacent atom
     */
    const int& getNatom(const int i) const
    {
        return adj_info.natom[i];
    }

    /**
     * @brief Get the position of an adjacent atom.
     * @param i index of adjacent atom
     * @return position of adjacent atom
     */
    const ModuleBase::Vector3<double>& getAdjacentTau(const int i) const
    {
        return adj_info.adjacent_tau[i];
    }

    /**
     * @brief Get the box indices of an adjacent atom.
     * @param i index of adjacent atom
     * @return box indices of adjacent atom
     */
    const ModuleBase::Vector3<int>& getBox(const int i) const
    {
        return adj_info.box[i];
    }

  private:
    mutable AdjacentAtomInfo adj_info; ///< adjacent atom information
    bool test_deconstructor;           ///< test deconstructor flag
};
#endif
