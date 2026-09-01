#include "source_pw/module_pwdft/onsite_proj.h"
#include "source_base/kernels/math_kernel_op.h"

template<typename T, typename Device>
projectors::OnsiteProjector<T, Device>* projectors::OnsiteProjector<T, Device>::get_instance()
{
    static projectors::OnsiteProjector<T, Device> instance;
    return &instance;
}

template<typename T, typename Device>
projectors::OnsiteProjector<T, Device>::~OnsiteProjector()
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
projectors::OnsiteProjector<double, base_device::DEVICE_CPU>*
projectors::OnsiteProjector<double, base_device::DEVICE_CPU>::get_instance();

template
projectors::OnsiteProjector<double, base_device::DEVICE_CPU>::~OnsiteProjector();

#if ((defined __CUDA) || (defined __ROCM))
template
projectors::OnsiteProjector<double, base_device::DEVICE_GPU>*
projectors::OnsiteProjector<double, base_device::DEVICE_GPU>::get_instance();

template
projectors::OnsiteProjector<double, base_device::DEVICE_GPU>::~OnsiteProjector();
#endif
