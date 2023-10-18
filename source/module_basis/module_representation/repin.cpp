#include "repin.h"

// Explicit instantiation
template class RepIn<std::complex<double>, psi::DEVICE_CPU>;
template class RepIn<std::complex<float>, psi::DEVICE_CPU>;
// gamma point support
template class RepIn<double, psi::DEVICE_CPU>;
template class RepIn<float, psi::DEVICE_CPU>;
// gpu support
#if ((defined __CUDA) || (defined __ROCM))
template class RepIn<std::complex<double>, psi::DEVICE_GPU>;
template class RepIn<std::complex<float>, psi::DEVICE_GPU>;
// gamma point support
template class RepIn<double, psi::DEVICE_GPU>;
template class RepIn<float, psi::DEVICE_GPU>;
#endif