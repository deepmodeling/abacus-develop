#include "source_base/vector3.h"
#include "source_cell/unitcell.h"
// mock of TD_info
class TD_info
{
  public:
    TD_info()
    {
    }
    ~TD_info()
    {
    }
    const UnitCell* get_ucell()
    {
        return nullptr;
    }
    static ModuleBase::Vector3<double> A_prop_ha;
    static TD_info* td_vel_op;
};
TD_info td_info;
TD_info* TD_info::td_vel_op = &td_info;
ModuleBase::Vector3<double> TD_info::A_prop_ha(0.0, 0.0, 0.0);
