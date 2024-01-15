#include "spin_constrain.h"

template <typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::set_adjs_all(const UnitCell& ucell, const Grid_Driver& GridD)
{
    std::cout << "set_adjs_all" << std::endl;
}

template class SpinConstrain<std::complex<double>, psi::DEVICE_CPU>;
template class SpinConstrain<double, psi::DEVICE_CPU>;