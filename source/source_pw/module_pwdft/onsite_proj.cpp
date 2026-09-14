#include "source_pw/module_pwdft/onsite_proj.h"
#include "source_base/kernels/math_kernel_op.h"

template<typename T, typename Device>
projectors::OnsiteProjector<T, Device>* projectors::OnsiteProjector<T, Device>::get_instance()
{
    static projectors::OnsiteProjector<T, Device> instance;
    return &instance;
}

namespace projectors {

template<typename T, typename Device>
OnsiteProjector<T, Device>::~OnsiteProjector()
{
    //delete[] becp;
    delete fs_tools;
    delmem_complex_op()(this->tab_atomic_);
    if(this->device == base_device::GpuDevice)
    {
        delmem_complex_h_op()(this->h_becp);
    }
    delmem_complex_op()(this->becp);

}

// explicit method instantiation
template
OnsiteProjector<double, base_device::DEVICE_CPU>*
OnsiteProjector<double, base_device::DEVICE_CPU>::get_instance();

template
OnsiteProjector<double, base_device::DEVICE_CPU>::~OnsiteProjector();

#if ((defined __CUDA) || (defined __ROCM))
template
OnsiteProjector<double, base_device::DEVICE_GPU>*
OnsiteProjector<double, base_device::DEVICE_GPU>::get_instance();

template
OnsiteProjector<double, base_device::DEVICE_GPU>::~OnsiteProjector();
#endif

} // namespace projectors
