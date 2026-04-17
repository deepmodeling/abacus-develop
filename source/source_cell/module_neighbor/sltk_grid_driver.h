#ifndef GRID_DRIVER_H
#define GRID_DRIVER_H

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/vector3.h"
#include "source_cell/unitcell.h"
#include "source_pw/module_pwdft/structure_factor.h"
#include "sltk_atom.h"
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
    Grid_Driver(){ test_deconstructor = false; };
    Grid_Driver(const int& test_d_in, const int& test_grid_in);

    ~Grid_Driver();

    Grid_Driver& operator=(Grid_Driver&&) = default;
    // [FSSH修改] 原代码: 无拷贝赋值运算符（因显式定义了移动赋值，拷贝赋值被隐式删除）
    // 修改原因: FSSH借用构造函数中需要拷贝Grid_Driver以保留KS对象的完整性，
    //   Grid_Driver及其父类Grid的所有成员均为可拷贝类型，default即可正确工作
    Grid_Driver& operator=(const Grid_Driver&) = default;

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
                   const int ntype,
                   const int nnumber,
                   AdjacentAtomInfo* adjs = nullptr) const;

    // cartesian_posi and ucell is deprecated 20241204 zhanghaochong
    // this interface is deprecated, please use Find_atom above
    void Find_atom(const UnitCell& ucell,
                   const ModuleBase::Vector3<double>& cartesian_posi,
                   const int& ntype,
                   const int& nnumber,
                   AdjacentAtomInfo* adjs = nullptr) const;
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
    bool test_deconstructor;
};
#endif
