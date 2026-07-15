#ifndef MD_STATE_VIEW_H
#define MD_STATE_VIEW_H

#include "source_base/vector3.h"
#include "source_cell/module_neighlist/neighbor_types.h"

#include <vector>

class MdCell;
class UnitCell;

struct MdAtomRef
{
    ModuleBase::Vector3<double>* cart = 0;
    ModuleBase::Vector3<double>* frac = 0;
    ModuleBase::Vector3<double>* vel = 0;
    ModuleBase::Vector3<int>* mbl = 0;
    ModuleBase::Vector3<double>* force = 0;
    double mass = 0.0;
    int type = 0;
    int type_index = 0;
    ModuleNeighList::GlobalAtomId global_id = -1;
};

class MdStateView
{
  public:
    static MdStateView from_unitcell(UnitCell& ucell,
                                     ModuleBase::Vector3<double>* force = 0,
                                     ModuleBase::Vector3<double>* vel = 0,
                                     ModuleBase::Vector3<int>* mbl = 0,
                                     double* mass = 0);
    static MdStateView from_mdcell(MdCell& mdcell);

    int size() const;
    bool uses_distributed_storage() const;

    ModuleBase::Vector3<double>& cart(const int i) const;
    ModuleBase::Vector3<double>& frac(const int i) const;
    ModuleBase::Vector3<double>& vel(const int i) const;
    ModuleBase::Vector3<int>& mbl(const int i) const;
    ModuleBase::Vector3<double>& force(const int i) const;

    double mass(const int i) const;
    int type(const int i) const;
    int type_index(const int i) const;
    ModuleNeighList::GlobalAtomId global_id(const int i) const;

  private:
    std::vector<MdAtomRef> atoms_;
    bool distributed_ = false;
};

#endif
