#ifndef VDW_PARAMETERS_H
#define VDW_PARAMETERS_H

#include "module_base/vector3.h"

#include <string>

namespace vdw
{

class VdwParameters
{
  public:
    VdwParameters() = default;
    virtual ~VdwParameters() = default;

    inline const std::string &model() const { return model_; }
    inline const ModuleBase::Vector3<int> &period() const { return period_; };

  protected:
    std::string model_;
    ModuleBase::Vector3<int> period_;
};

} // namespace vdw

#endif // VDW_PARAMETERS_H
