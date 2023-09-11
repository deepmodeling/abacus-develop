#include "veff_pw.h"

#include "module_base/timer.h"
#include "module_base/tool_quit.h"
#include "module_psi/kernels/device.h"

using hamilt::Veff;
using hamilt::OperatorPW;

template <typename FPTYPE, typename Device>
Veff<OperatorPW<FPTYPE, Device>>::Veff(const int* isk_in,
                                       const FPTYPE* veff_in,
                                       const int veff_row,
                                       const int veff_col,
                                       const ModulePW::PW_Basis_K* wfcpw_in)
{
    this->classname = "Veff";
    this->cal_type = pw_veff;
    this->isk = isk_in;
    this->veff = veff_in;
    //note: "veff = nullptr" means that this core does not treat potential but still treats wf. 
    this->veff_row = veff_row;
    this->veff_col = veff_col;
    this->wfcpw = wfcpw_in;
    resmem_complex_op()(this->ctx, this->porter, this->wfcpw->nmaxgr, "Veff<PW>::porter");
    resmem_complex_op()(this->ctx, this->porter1, this->wfcpw->nmaxgr, "Veff<PW>::porter1");
    if (this->isk == nullptr || this->wfcpw == nullptr) {
        ModuleBase::WARNING_QUIT("VeffPW", "Constuctor of Operator::VeffPW is failed, please check your code!");
    }
}

template<typename FPTYPE, typename Device>
Veff<OperatorPW<FPTYPE, Device>>::~Veff()
{
    delmem_complex_op()(this->ctx, this->porter);
    delmem_complex_op()(this->ctx, this->porter1);
}

template<typename FPTYPE, typename Device>
void Veff<OperatorPW<FPTYPE, Device>>::act(
    const int nbands,
    const int nbasis,
    const int npol,
    const std::complex<FPTYPE>* tmpsi_in,
    std::complex<FPTYPE>* tmhpsi,
    const int ngk_ik)const
{
    ModuleBase::timer::tick("Operator", "VeffPW");

    int max_npw = nbasis / npol;
    const int current_spin = this->isk[this->ik];
    
    // std::complex<FPTYPE> *porter = new std::complex<FPTYPE>[wfcpw->nmaxgr];
    for (int ib = 0; ib < nbands; ib += npol)
    {
        if (npol == 1)
        {
            // wfcpw->recip2real(tmpsi_in, porter, this->ik);
            wfcpw->recip_to_real(this->ctx, tmpsi_in, this->porter, this->ik);
            // NOTICE: when MPI threads are larger than number of Z grids
            // veff would contain nothing, and nothing should be done in real space
            // but the 3DFFT can not be skipped, it will cause hanging
            if(this->veff_col != 0)
            {
                veff_op()(this->ctx, this->veff_col, this->porter, this->veff + current_spin * this->veff_col);
                // const FPTYPE* current_veff = &(this->veff[0](current_spin, 0));
                // for (int ir = 0; ir < this->veff->nc; ++ir)
                // {
                //     porter[ir] *= current_veff[ir];
                // }
            }
            // wfcpw->real2recip(porter, tmhpsi, this->ik, true);
            wfcpw->real_to_recip(this->ctx, this->porter, tmhpsi, this->ik, true);
        }
        else
        {
            // std::complex<FPTYPE> *porter1 = new std::complex<FPTYPE>[wfcpw->nmaxgr];
            // fft to real space and doing things.
            wfcpw->recip_to_real(this->ctx, tmpsi_in, this->porter, this->ik);
            wfcpw->recip_to_real(this->ctx, tmpsi_in + max_npw, this->porter1, this->ik);
            if(this->veff_col != 0)
            {
                /// denghui added at 20221109
                const FPTYPE* current_veff[4];
                for(int is = 0; is < 4; is++) {
                    current_veff[is] = this->veff + is * this->veff_col ; // for CPU device
                }
                veff_op()(this->ctx, this->veff_col, this->porter, this->porter1, current_veff);
                // std::complex<FPTYPE> sup, sdown;
                // for (int ir = 0; ir < this->veff_col; ir++) {
                //     sup = this->porter[ir] * (current_veff[0][ir] + current_veff[3][ir])
                //         + this->porter1[ir]
                //                 * (current_veff[1][ir]
                //                 - std::complex<FPTYPE>(0.0, 1.0) * current_veff[2][ir]);
                //     sdown = this->porter1[ir] * (current_veff[0][ir] - current_veff[3][ir])
                //             + this->porter[ir]
                //                 * (current_veff[1][ir]
                //                     + std::complex<FPTYPE>(0.0, 1.0) * current_veff[2][ir]);
                //     this->porter[ir] = sup;
                //     this->porter1[ir] = sdown;
                // }
            }
            // (3) fft back to G space.
            wfcpw->real_to_recip(this->ctx, this->porter, tmhpsi, this->ik, true);
            wfcpw->real_to_recip(this->ctx, this->porter1, tmhpsi + max_npw, this->ik, true);
        }
        tmhpsi += max_npw * npol;
        tmpsi_in += max_npw * npol;
    }
    ModuleBase::timer::tick("Operator", "VeffPW");
}

template<typename FPTYPE, typename Device>
template<typename T_in, typename Device_in>
hamilt::Veff<OperatorPW<FPTYPE, Device>>::Veff(const Veff<OperatorPW<T_in, Device_in>> *veff) {
    this->classname = "Veff";
    this->cal_type = pw_veff;
    this->ik = veff->get_ik();
    this->isk = veff->get_isk();
    this->veff_col = veff->get_veff_col();
    this->veff_row = veff->get_veff_row();
    this->wfcpw = veff->get_wfcpw();
    resmem_complex_op()(this->ctx, this->porter, this->wfcpw->nmaxgr);
    resmem_complex_op()(this->ctx, this->porter1, this->wfcpw->nmaxgr);
    this->veff = veff->get_veff();
    if (this->isk == nullptr || this->veff == nullptr || this->wfcpw == nullptr) {
        ModuleBase::WARNING_QUIT("VeffPW", "Constuctor of Operator::VeffPW is failed, please check your code!");
    }
}

namespace hamilt {
template class Veff<OperatorPW<float, psi::DEVICE_CPU>>;
template class Veff<OperatorPW<double, psi::DEVICE_CPU>>;
// template Veff<OperatorPW<double, psi::DEVICE_CPU>>::Veff(const Veff<OperatorPW<double, psi::DEVICE_CPU>> *veff);
#if ((defined __CUDA) || (defined __ROCM))
template class Veff<OperatorPW<float, psi::DEVICE_GPU>>;
template class Veff<OperatorPW<double, psi::DEVICE_GPU>>;
// template Veff<OperatorPW<double, psi::DEVICE_GPU>>::Veff(const Veff<OperatorPW<double, psi::DEVICE_GPU>> *veff);
#endif
} // namespace hamilt