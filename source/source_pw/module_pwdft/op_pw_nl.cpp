#include "op_pw_nl.h"

#include "source_base/timer.h"
#include "source_base/tool_quit.h"

namespace hamilt
{

template <typename T, typename Device>
Nonlocal<OperatorPW<T, Device>>::Nonlocal(const int* isk_in,
                                          const pseudopot_cell_vnl* ppcell_in,
                                          const UnitCell* ucell_in,
                                          const ModulePW::PW_Basis_K* wfc_basis)
{
    if (isk_in == nullptr || ppcell_in == nullptr || ucell_in == nullptr || wfc_basis == nullptr)
    {
        ModuleBase::WARNING_QUIT("NonlocalPW", "Constuctor of Operator::NonlocalPW is failed, please check your code!");
    }
    this->classname = "Nonlocal";
    this->cal_type = calculation_type::pw_nonlocal;
    this->wfcpw = wfc_basis;
    this->isk = isk_in;
    this->ppcell = ppcell_in;
    this->ucell = ucell_in;
    this->vkb = this->ppcell->template get_vkb_data<Real>();
}

template <typename T, typename Device>
Nonlocal<OperatorPW<T, Device>>::~Nonlocal()
{
}

template <typename T, typename Device>
void Nonlocal<OperatorPW<T, Device>>::init(const int ik_in)
{
    ModuleBase::timer::start("Nonlocal", "init");
    this->ik = ik_in;
    // Calculate nonlocal pseudopotential vkb
    if (this->ppcell->nkb > 0) // xiaohui add 2013-09-02. Attention...
    {
        this->ppcell->getvnl(this->ctx, *this->ucell, this->ik, ModuleBase::Vector3<double>(0.0, 0.0, 0.0), this->vkb);
    }

    if (this->next_op != nullptr)
    {
        this->next_op->init(ik_in);
    }

    ModuleBase::timer::end("Nonlocal", "init");
}

template <typename T, typename Device>
void Nonlocal<OperatorPW<T, Device>>::act(const int nbands,
                                          const int nbasis,
                                          const int npol,
                                          const T* tmpsi_in,
                                          T* tmhpsi,
                                          const int ngk_ik,
                                          const bool is_first_node) const
{
    this->workspace_.apply(*this->ucell,
                           *this->ppcell,
                           *this->wfcpw,
                           this->isk[this->ik],
                           nbands,
                           nbasis,
                           npol,
                           ngk_ik,
                           is_first_node,
                           this->vkb,
                           tmpsi_in,
                           tmhpsi);
}

template <typename T, typename Device>
template <typename T_in, typename Device_in>
hamilt::Nonlocal<OperatorPW<T, Device>>::Nonlocal(const Nonlocal<OperatorPW<T_in, Device_in>>* nonlocal)
{
    this->classname = "Nonlocal";
    this->cal_type = calculation_type::pw_nonlocal;
    this->ik = nonlocal->get_ik();
    this->isk = nonlocal->get_isk();
    this->ppcell = nonlocal->get_ppcell();
    this->ucell = nonlocal->get_ucell();
    this->wfcpw = nonlocal->get_wfcpw();
    this->vkb = this->ppcell->template get_vkb_data<Real>();
    if (this->isk == nullptr || this->ppcell == nullptr || this->ucell == nullptr || this->wfcpw == nullptr)
    {
        ModuleBase::WARNING_QUIT("NonlocalPW", "Constuctor of Operator::NonlocalPW is failed, please check your code!");
    }
}

template class Nonlocal<OperatorPW<std::complex<float>, base_device::DEVICE_CPU>>;
template class Nonlocal<OperatorPW<std::complex<double>, base_device::DEVICE_CPU>>;
// template Nonlocal<OperatorPW<std::complex<double>, base_device::DEVICE_CPU>>::Nonlocal(const
// Nonlocal<OperatorPW<std::complex<double>, base_device::DEVICE_CPU>> *nonlocal);
#if ((defined __CUDA) || (defined __ROCM))
template class Nonlocal<OperatorPW<std::complex<float>, base_device::DEVICE_GPU>>;
template class Nonlocal<OperatorPW<std::complex<double>, base_device::DEVICE_GPU>>;
// template Nonlocal<OperatorPW<std::complex<double>, base_device::DEVICE_CPU>>::Nonlocal(const
// Nonlocal<OperatorPW<std::complex<double>, base_device::DEVICE_GPU>> *nonlocal); template
// Nonlocal<OperatorPW<std::complex<double>, base_device::DEVICE_GPU>>::Nonlocal(const
// Nonlocal<OperatorPW<std::complex<double>, base_device::DEVICE_CPU>> *nonlocal); template
// Nonlocal<OperatorPW<std::complex<double>, base_device::DEVICE_GPU>>::Nonlocal(const
// Nonlocal<OperatorPW<std::complex<double>, base_device::DEVICE_GPU>> *nonlocal);
#endif
} // namespace hamilt
