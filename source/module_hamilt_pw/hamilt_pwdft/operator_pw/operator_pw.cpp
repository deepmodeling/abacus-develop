#include "module_base/timer.h"
#include "module_hamilt_general/operator.h"
#include "module_hamilt_pw/hamilt_pwdft/operator_pw/operator_pw.h"

using namespace hamilt;

template<typename T, typename Device>
OperatorPW<T, Device>::~OperatorPW(){};

namespace hamilt {
template class OperatorPW<std::complex<float>, psi::DEVICE_CPU>;
template class OperatorPW<std::complex<double>, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class OperatorPW<std::complex<float>, psi::DEVICE_GPU>;
template class OperatorPW<std::complex<double>, psi::DEVICE_GPU>;
#endif
}