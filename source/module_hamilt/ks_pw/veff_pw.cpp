#include "veff_pw.h"

#include "module_base/timer.h"
#include "module_base/tool_quit.h"

using hamilt::Veff;
using hamilt::OperatorPW;

template<typename FPTYPE, typename Device>
Veff<OperatorPW<FPTYPE, Device>>::Veff(
    const int* isk_in,
    const ModuleBase::matrix* veff_in,
    ModulePW::PW_Basis_K* wfcpw_in)
{
    this->cal_type = pw_veff;
    this->isk = isk_in;
    // this->veff = veff_in;
    // TODO: add an GPU veff array
    this->veff = veff_in[0].c;
    this->veff_col = veff_in[0].nc;
    this->wfcpw = wfcpw_in;
    resize_memory_op()(this->ctx, this->porter, this->wfcpw->nmaxgr);
    if (this->npol != 1) {
        resize_memory_op()(this->ctx, this->porter1, this->wfcpw->nmaxgr);
    }
    if (this->isk == nullptr || this->veff == nullptr || this->wfcpw == nullptr) {
        ModuleBase::WARNING_QUIT("VeffPW", "Constuctor of Operator::VeffPW is failed, please check your code!");
    }
}

template<typename FPTYPE, typename Device>
Veff<OperatorPW<FPTYPE, Device>>::~Veff()
{
    delete_memory_op()(this->ctx, this->porter);
    if (this->npol != 1) {
        delete_memory_op()(this->ctx, this->porter1);
    }
}

template<typename FPTYPE, typename Device>
void Veff<OperatorPW<FPTYPE, Device>>::act(
    const psi::Psi<std::complex<FPTYPE>, Device> *psi_in, 
    const int n_npwx, 
    const std::complex<FPTYPE>* tmpsi_in, 
    std::complex<FPTYPE>* tmhpsi
)const  
{
    ModuleBase::timer::tick("Operator", "VeffPW");

    this->max_npw = psi_in->get_nbasis() / psi_in->npol;
    const int current_spin = this->isk[this->ik];
    this->npol = psi_in->npol;
    
    // std::complex<FPTYPE> *porter = new std::complex<FPTYPE>[wfcpw->nmaxgr];
    for (int ib = 0; ib < n_npwx; ib += this->npol)
    {
        if (this->npol == 1)
        {
            // wfcpw->recip2real(tmpsi_in, porter, this->ik);
            wfcpw->recip_to_real(this->ctx, tmpsi_in, this->porter, this->ik);
            // NOTICE: when MPI threads are larger than number of Z grids
            // veff would contain nothing, and nothing should be done in real space
            // but the 3DFFT can not be skipped, it will cause hanging
            if(this->veff_col != 0)
            {
                // const FPTYPE* current_veff = &(this->veff[0](current_spin, 0));
                // for (int ir = 0; ir < this->veff->nc; ++ir)
                // {
                //     porter[ir] *= current_veff[ir];
                // }
                veff_op()(this->ctx, this->veff_col, this->porter, this->veff + current_spin * this->veff_col);
            }
            // wfcpw->real2recip(porter, tmhpsi, this->ik, true);
            wfcpw->real_to_recip(this->ctx, this->porter, tmhpsi, this->ik, true);
        }
        else
        {
            // std::complex<FPTYPE> *porter1 = new std::complex<FPTYPE>[wfcpw->nmaxgr];
            // fft to real space and doing things.
            wfcpw->recip2real(tmpsi_in, this->porter, this->ik);
            wfcpw->recip2real(tmpsi_in + this->max_npw, this->porter1, this->ik);
            std::complex<FPTYPE> sup, sdown;
            if(this->veff_col != 0)
            {
                const FPTYPE* current_veff[4];
                for(int is=0;is<4;is++)
                {
                    current_veff[is] = this->veff + is * this->veff_col;
                }
                for (int ir = 0; ir < this->veff_col; ir++)
                {
                    sup = this->porter[ir] * (current_veff[0][ir] + current_veff[3][ir])
                        + this->porter1[ir]
                                * (current_veff[1][ir]
                                - std::complex<FPTYPE>(0.0, 1.0) * current_veff[2][ir]);
                    sdown = this->porter1[ir] * (current_veff[0][ir] - current_veff[3][ir])
                            + this->porter[ir]
                                * (current_veff[1][ir]
                                    + std::complex<FPTYPE>(0.0, 1.0) * current_veff[2][ir]);
                    this->porter[ir] = sup;
                    this->porter1[ir] = sdown;
                }
            }
            // (3) fft back to G space.
            wfcpw->real2recip(this->porter, tmhpsi, this->ik, true);
            wfcpw->real2recip(this->porter1, tmhpsi + this->max_npw, this->ik, true);
        }
        tmhpsi += this->max_npw * this->npol;
        tmpsi_in += this->max_npw * this->npol;
    }
    ModuleBase::timer::tick("Operator", "VeffPW");
}

namespace hamilt{
template class Veff<OperatorPW<double, psi::DEVICE_CPU>>;
#if ((defined __CUDA) || (defined __ROCM))
template class Veff<OperatorPW<double, psi::DEVICE_GPU>>;
#endif
} // namespace hamilt