#include "spin_constrain.h"

template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::cal_mw_from_lambda()
{
    std::cout << "hello world" << std::endl;
}

template class SpinConstrain<double, psi::DEVICE_CPU>;