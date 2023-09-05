#include "module_base/timer.h"
#include "module_hamilt_general/operator.h"
#include "module_hamilt_pw/hamilt_pwdft/operator_pw/operator_pw.h"

using namespace hamilt;

template<typename FPTYPE, typename Device>
OperatorPW<FPTYPE, Device>::~OperatorPW(){};

namespace hamilt {
template class OperatorPW<float, psi::DEVICE_CPU>;
template class OperatorPW<double, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class OperatorPW<float, psi::DEVICE_GPU>;
template class OperatorPW<double, psi::DEVICE_GPU>;
#endif
}