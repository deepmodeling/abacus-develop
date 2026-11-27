#include "source_lcao/module_rt/td_info.h"

// These DeePKS tests do not create a TD context.
const ModuleBase::Vector3<double>& TD_info::endpoint_A() const
{
    static const ModuleBase::Vector3<double> value(0.0, 0.0, 0.0);
    return value;
}

int elecstate::TDFieldManager::gauge() const
{
    return gauge_;
}
