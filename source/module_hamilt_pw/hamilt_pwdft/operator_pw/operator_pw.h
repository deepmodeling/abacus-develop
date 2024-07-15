#ifndef OPERATORPW_H
#define OPERATORPW_H
#include"module_hamilt_general/operator.h"

namespace hamilt {
template <typename T, typename Device = base_device::DEVICE_CPU>
class OperatorPW : public Operator<T, Device>
{
    public:
        virtual ~OperatorPW();
        std::string classname = "";
};

}//end namespace hamilt

#endif