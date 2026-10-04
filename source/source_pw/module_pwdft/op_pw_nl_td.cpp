#include "op_pw_nl_td.h"

#include "source_base/timer.h"
#include "source_base/tool_quit.h"

namespace hamilt
{

template <typename T, typename Device>
TDNonlocalPW<T, Device>::TDNonlocalPW(const int* isk_in,
                                      const pseudopot_cell_vnl* ppcell_in,
                                      const UnitCell* ucell_in,
                                      const ModulePW::PW_Basis_K* wfc_basis)
{
    if (isk_in == nullptr || ppcell_in == nullptr || ucell_in == nullptr || wfc_basis == nullptr)
    {
        ModuleBase::WARNING_QUIT("TDNonlocalPW", "Constructor failed, null pointers detected!");
    }

    this->classname = "TDNonlocalPW";
    this->cal_type = calculation_type::pw_nonlocal;
    this->isk = isk_in;
    this->ppcell = ppcell_in;
    this->ucell = ucell_in;
    this->wfcpw = wfc_basis;
}

template <typename T, typename Device>
TDNonlocalPW<T, Device>::~TDNonlocalPW()
{
    delmem_complex_op()(this->vkb_td);
}

template <typename T, typename Device>
void TDNonlocalPW<T, Device>::init(const int ik_in)
{
    ModuleBase::timer::start("TDNonlocalPW", "init");
    this->ik = ik_in;

    // Refresh the shifted projectors when nonlocal projectors are present.
    if (this->ppcell->nkb > 0 && this->wfcpw->npwk[ik_in] > 0)
    {
        // Allocate the time-dependent projector cache on first use.
        if (this->vkb_td == nullptr)
        {
            resmem_complex_op()(this->vkb_td, this->ppcell->nkb * this->wfcpw->npwk_max, "TDNL::vkb_td");
        }

        // Hartree-unit vector potential.
        ModuleBase::Vector3<double> A_au = A_ha_;

        // Generate projectors shifted by the current vector potential.
        this->ppcell->getvnl(this->ctx, *this->ucell, this->ik, A_au, this->vkb_td);
    }

    if (this->next_op != nullptr)
    {
        this->next_op->init(ik_in);
    }
    ModuleBase::timer::end("TDNonlocalPW", "init");
}

template <typename T, typename Device>
void TDNonlocalPW<T, Device>::act(const int nbands,
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
                           this->vkb_td,
                           tmpsi_in,
                           tmhpsi);
}

// Explicit CPU template instantiations.
template class TDNonlocalPW<std::complex<float>, base_device::DEVICE_CPU>;
template class TDNonlocalPW<std::complex<double>, base_device::DEVICE_CPU>;

// ================= GPU explicit instantiations =================
#if ((defined __CUDA) || (defined __ROCM))
template class TDNonlocalPW<std::complex<float>, base_device::DEVICE_GPU>;
template class TDNonlocalPW<std::complex<double>, base_device::DEVICE_GPU>;
#endif
// ===============================================================

} // namespace hamilt
