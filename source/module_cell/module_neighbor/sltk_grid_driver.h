#ifndef GRID_DRIVER_H
#define GRID_DRIVER_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/vector3.h"
#include "module_cell/unitcell.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "sltk_atom.h"
#include "sltk_atom_input.h"
#include "sltk_grid.h"

#include <memory>
#include <stdexcept>
#include <tuple>

//==========================================================
// Struct of array for packing the Adjacent atom information
//==========================================================
class AdjacentAtomInfo
{
  public:
    AdjacentAtomInfo() : adj_num(0)
    {
    }
    int adj_num;
    std::vector<int> ntype;
    std::vector<int> natom;
    std::vector<ModuleBase::Vector3<double>> adjacent_tau;
    std::vector<ModuleBase::Vector3<int>> box;
    void clear()
    {
        adj_num = 0;
        ntype.clear();
        natom.clear();
        adjacent_tau.clear();
        box.clear();
    }
};

void filter_adjs(const std::vector<bool>& is_adj, AdjacentAtomInfo& adjs);

class Grid_Driver : public Grid
{
  public:
    //==========================================================
    // THE INTERFACE WITH USER :
    // MEMBRE FUNCTIONS :
    // NAME : Find_atom (input cartesian position,find the
    //		adjacent of this atom,and store the information
    //		in 'adj_num','ntype','natom'
    //==========================================================
    Grid_Driver() : test_deconstructor(0){};
    Grid_Driver(const int& test_d_in, const int& test_grid_in);

    ~Grid_Driver();

    //==========================================================
    // EXPLAIN FOR default parameter `adjs = nullptr`
    //
    // This design make Grid_Driver compatible with multi-thread usage
    // 1. Find_atom store results in Grid_Driver::adj_info
    //     by default.
    // 2. And store results into parameter adjs when adjs is
    //     NOT NULL
    //==========================================================
    void Find_atom(const UnitCell& ucell,
                   const ModuleBase::Vector3<double>& cartesian_posi,
                   const int& ntype,
                   const int& nnumber,
                   AdjacentAtomInfo* adjs = nullptr);

    //==========================================================
    // EXPLAIN : The adjacent information for the input
    // cartesian_pos
    // MEMBER VARIABLES :
    // NAME : getAdjacentNum
    // NAME : getNtype
    // NAME : getNatom
    // NAME : getAdjaentTau
    //==========================================================
    const int& getAdjacentNum() const
    {
        return adj_info.adj_num;
    }
    const int& getType(const int i) const
    {
        return adj_info.ntype[i];
    }
    const int& getNatom(const int i) const
    {
        return adj_info.natom[i];
    }
    const ModuleBase::Vector3<double>& getAdjacentTau(const int i) const
    {
        return adj_info.adjacent_tau[i];
    }
    const ModuleBase::Vector3<int>& getBox(const int i) const
    {
        return adj_info.box[i];
    }

  private:
    mutable AdjacentAtomInfo adj_info;

    const int test_deconstructor; // caoyu reconst 2021-05-24

    //==========================================================
    // MEMBER FUNCTIONS :
    // NAME : Calculate_adjacent_site
    //==========================================================
    ModuleBase::Vector3<double> Calculate_adjacent_site(const double x,
                                                        const double y,
                                                        const double z,
                                                        const double& box11,
                                                        const double& box12,
                                                        const double& box13,
                                                        const double& box21,
                                                        const double& box22,
                                                        const double& box23,
                                                        const double& box31,
                                                        const double& box32,
                                                        const double& box33,
                                                        const short box_x, // three dimensions of the target box
                                                        const short box_y,
                                                        const short box_z) const;
};

namespace GlobalC
{
extern Grid_Driver GridD;
}
#endif
