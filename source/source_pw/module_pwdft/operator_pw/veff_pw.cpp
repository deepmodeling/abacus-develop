#include "veff_pw.h"

#include "source_base/timer.h"
#include "source_base/tool_quit.h"

namespace hamilt {

template<typename T, typename Device>
Veff<OperatorPW<T, Device>>::Veff(const int* isk_in,
                                       const Real* veff_in,
                                       const int veff_row,
                                       const int veff_col,
                                       const ModulePW::PW_Basis_K* wfcpw_in)
{
    if (isk_in == nullptr || wfcpw_in == nullptr) 
    {
        ModuleBase::WARNING_QUIT("VeffPW", "Constuctor of Operator::VeffPW is failed, please check your code!");
    }

    this->classname = "Veff";
    this->cal_type = calculation_type::pw_veff;
    this->isk = isk_in;
    this->veff = veff_in;
    //note: "veff = nullptr" means that this core does not treat potential but still treats wf. 
    this->veff_row = veff_row;
    this->veff_col = veff_col;
    this->wfcpw = wfcpw_in;
    resmem_complex_op()(this->porter, this->wfcpw->nmaxgr, "Veff<PW>::porter");
    resmem_complex_op()(this->porter1, this->wfcpw->nmaxgr, "Veff<PW>::porter1");

}

template<typename T, typename Device>
Veff<OperatorPW<T, Device>>::~Veff()
{
    delmem_complex_op()(this->porter);
    delmem_complex_op()(this->porter1);
}

template<typename T, typename Device>
void Veff<OperatorPW<T, Device>>::act(
    const int nbands,
    const int nbasis,
    const int npol,
    const T* tmpsi_in,
    T* tmhpsi,
    const int ngk_ik,
    const bool is_first_node)const
{
    ModuleBase::timer::tick("Operator", "veff_pw");
    if(is_first_node)
    {
        setmem_complex_op()(tmhpsi, 0, nbasis*nbands/npol);
    }
    int max_npw = nbasis / npol;
    const int current_spin = this->isk[this->ik];
    const int psi_offset= max_npw * npol;
#ifdef __DSP
    if (npol == 1)
    {
        ModulePW::FFT_Guard guard(wfcpw->fft_bundle);
        for (int ib = 0; ib < nbands; ib += npol)
        {
            wfcpw->convolution(this->ctx,
                               this->ik,
                               this->veff_col,
                               tmpsi_in,
                               this->veff + current_spin * this->veff_col,
                               tmhpsi,
                               true);
            tmhpsi   += psi_offset;
            tmpsi_in += psi_offset;
        }
    }else if (npol == 2)
    {
        const Real* current_veff[4]={nullptr};
        for (int is = 0; is < 4; is++)
        {
            current_veff[is] = this->veff + is * this->veff_col;
        }
        for (int ib = 0; ib < nbands; ib += npol)
        {
            wfcpw->recip_to_real<T, Device>(tmpsi_in, this->porter, this->ik);
            wfcpw->recip_to_real<T, Device>(tmpsi_in + max_npw, this->porter1, this->ik);
            veff_op()(this->ctx, this->veff_col, this->porter, this->porter1, current_veff);
            wfcpw->real_to_recip<T, Device>(this->porter, tmhpsi, this->ik, true);
            wfcpw->real_to_recip<T, Device>(this->porter1, tmhpsi + max_npw, this->ik, true);
            tmhpsi   += psi_offset;
            tmpsi_in += psi_offset;
        }
    }else{
        ModuleBase::WARNING_QUIT("VeffPW", "npol should be 1 or 2 or veff_col equal to 0\n");
    }
#else
    if (npol == 1)
    {
        for (int ib = 0; ib < nbands; ib += npol)
        {
            wfcpw->recip_to_real<T, Device>(tmpsi_in, this->porter, this->ik);
            // NOTICE: when MPI threads are larger than the number of Z grids
            // veff would contain nothing, and nothing should be done in real space
            // but the 3DFFT can not be skipped, it will cause hanging
            veff_op()(this->ctx, this->veff_col, this->porter, this->veff + current_spin * this->veff_col);
            wfcpw->real_to_recip<T, Device>(this->porter, tmhpsi, this->ik, true);
            tmhpsi   += psi_offset;
            tmpsi_in += psi_offset;
        }
    }
    else if (npol == 2)
    {
        const Real* current_veff[4]={nullptr};
        for (int is = 0; is < 4; is++)
        {
            current_veff[is] = this->veff + is * this->veff_col;
        }
        for (int ib = 0; ib < nbands; ib += npol)
        {
            // FFT to real space and do things.
            wfcpw->recip_to_real<T, Device>(tmpsi_in, this->porter, this->ik);
            wfcpw->recip_to_real<T, Device>(tmpsi_in + max_npw, this->porter1, this->ik);
            veff_op()(this->ctx, this->veff_col, this->porter, this->porter1, current_veff);
            // FFT back to G space.
            wfcpw->real_to_recip<T, Device>(this->porter, tmhpsi, this->ik, true);
            wfcpw->real_to_recip<T, Device>(this->porter1, tmhpsi + max_npw, this->ik, true);
            tmhpsi   += psi_offset;
            tmpsi_in += psi_offset;
        }
    }else{
        ModuleBase::WARNING_QUIT("VeffPW", "npol should be 1 or 2 or veff_col equal to 0\n");
    }
#endif
    ModuleBase::timer::tick("Operator", "veff_pw");
}

template<typename T, typename Device>
template<typename T_in, typename Device_in>
hamilt::Veff<OperatorPW<T, Device>>::Veff(const Veff<OperatorPW<T_in, Device_in>> *veff) {
    this->classname = "Veff";
    this->cal_type = calculation_type::pw_veff;
    this->ik = veff->get_ik();
    this->isk = veff->get_isk();
    this->veff_col = veff->get_veff_col();
    this->veff_row = veff->get_veff_row();
    this->wfcpw = veff->get_wfcpw();
    resmem_complex_op()(this->porter, this->wfcpw->nmaxgr);
    resmem_complex_op()(this->porter1, this->wfcpw->nmaxgr);
    this->veff = veff->get_veff();
    if (this->isk == nullptr || this->veff == nullptr || this->wfcpw == nullptr) {
        ModuleBase::WARNING_QUIT("VeffPW", "Constuctor of Operator::VeffPW is failed, please check your code!");
    }
}

template class Veff<OperatorPW<std::complex<float>, base_device::DEVICE_CPU>>;
template class Veff<OperatorPW<std::complex<double>, base_device::DEVICE_CPU>>;
// template Veff<OperatorPW<std::complex<double>, base_device::DEVICE_CPU>>::Veff(const
// Veff<OperatorPW<std::complex<double>, base_device::DEVICE_CPU>> *veff);
#if ((defined __CUDA) || (defined __ROCM))
template class Veff<OperatorPW<std::complex<float>, base_device::DEVICE_GPU>>;
template class Veff<OperatorPW<std::complex<double>, base_device::DEVICE_GPU>>;
// template Veff<OperatorPW<std::complex<double>, base_device::DEVICE_GPU>>::Veff(const
// Veff<OperatorPW<std::complex<double>, base_device::DEVICE_GPU>> *veff);
#endif
} // namespace hamilt
