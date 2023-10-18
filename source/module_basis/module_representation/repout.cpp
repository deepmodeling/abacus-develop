#include "repout.h"

template<typename T, typename Device>
RepOut<T, Device>::RepOut()
{
}

template<typename T, typename Device>
RepOut<T, Device>::~RepOut()
{
}

// Explicit instantiation
template class RepOut<std::complex<double>, psi::DEVICE_CPU>;
template class RepOut<std::complex<float>, psi::DEVICE_CPU>;
// gamma point support
template class RepOut<double, psi::DEVICE_CPU>;
template class RepOut<float, psi::DEVICE_CPU>;
// gpu support
#if ((defined __CUDA) || (defined __ROCM))
template class RepOut<std::complex<double>, psi::DEVICE_GPU>;
template class RepOut<std::complex<float>, psi::DEVICE_GPU>;
// gamma point support
template class RepOut<double, psi::DEVICE_GPU>;
template class RepOut<float, psi::DEVICE_GPU>;
#endif