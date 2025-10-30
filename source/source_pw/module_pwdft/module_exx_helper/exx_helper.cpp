#include "exx_helper.h"

template <typename T, typename Device>
double Exx_Helper<T, Device>::cal_exx_energy(psi::Psi<T, Device> *psi_)
{
    return op_exx->cal_exx_energy(psi_);

}

template <typename T, typename Device>
bool Exx_Helper<T, Device>::exx_after_converge(int &iter, bool ene_conv)
{
    if (op_exx->first_iter)
    {
        op_exx->first_iter = false;
    }
    else if (!GlobalC::exx_info.info_global.separate_loop)
    {
        return true;
    }
    else if (PARAM.inp.exx_thr_type == "energy" && ene_conv)
    {
        return true;
    }
    else if (PARAM.inp.exx_thr_type == "density" && iter == 1)
    {
        return true;
    }
    else if (iter >= PARAM.inp.exx_hybrid_step)
    {
        GlobalV::ofs_running << " !!EXX IS NOT CONVERGED!!" << std::endl;
        std::cout << " !!EXX IS NOT CONVERGED!!" << std::endl;
        return true;
    }
    GlobalV::ofs_running << "Updating EXX and rerun SCF" << std::endl;
    iter = 0;
    return false;

}

template <typename T, typename Device>
void Exx_Helper<T, Device>::set_psi(psi::Psi<T, Device> *psi_)
{
    if (psi_ == nullptr)
        return;
    op_exx->set_psi(*psi_);
    if (PARAM.inp.exxace && GlobalC::exx_info.info_global.separate_loop)
    {
        op_exx->construct_ace();
    }
}

template class Exx_Helper<std::complex<float>, base_device::DEVICE_CPU>;
template class Exx_Helper<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Exx_Helper<std::complex<float>, base_device::DEVICE_GPU>;
template class Exx_Helper<std::complex<double>, base_device::DEVICE_GPU>;
#endif

#ifndef __EXX
#include "source_hamilt/module_xc/exx_info.h"
namespace GlobalC
{
    Exx_Info exx_info;
}
#endif