#include "veff_pw.h"

#include "module_base/timer.h"
#include "src_pw/global.h"
#include "module_base/tool_quit.h"

namespace hamilt
{

Veff<OperatorPW>::Veff(
    const int* isk_in,
    const ModuleBase::matrix* veff_in,
    ModulePW::PW_Basis_K* wfcpw_in
)
{
    this->cal_type = 11;
    this->isk = isk_in;
    this->veff = veff_in;
    this->wfcpw = wfcpw_in;
    if( this->isk == nullptr || this->veff == nullptr || this->wfcpw == nullptr)
    {
        ModuleBase::WARNING_QUIT("VeffPW", "Constuctor of Operator::VeffPW is failed, please check your code!");
    }
}

void Veff<OperatorPW>::act
(
    const psi::Psi<std::complex<double>> *psi_in, 
    const int n_npwx, 
    const std::complex<double>* tmpsi_in, 
    std::complex<double>* tmhpsi
)const  
{
    ModuleBase::timer::tick("Operator", "VeffPW");

    this->max_npw = psi_in->get_nbasis() / psi_in->npol;
    const int current_spin = this->isk[this->ik];
    this->npol = psi_in->npol;
    
    std::complex<double> *porter = new std::complex<double>[wfcpw->nmaxgr];
    for (int ib = 0; ib < n_npwx; ib += this->npol)
    {
        if (this->npol == 1)
        {
            wfcpw->recip2real(tmpsi_in, porter, ik);
            // NOTICE: when MPI threads are larger than number of Z grids
            // veff would contain nothing, and nothing should be done in real space
            // but the 3DFFT can not be skipped, it will cause hanging
            if(this->veff->nc != 0)
            {
                const double* current_veff = &(this->veff[0](current_spin, 0));
                for (int ir = 0; ir < this->veff->nc; ++ir)
                {
                    porter[ir] *= current_veff[ir];
                }
            }
            wfcpw->real2recip(porter, tmhpsi, ik, true);
        }
        else
        {
            std::complex<double> *porter1 = new std::complex<double>[wfcpw->nmaxgr];
            // fft to real space and doing things.
            wfcpw->recip2real(tmpsi_in, porter, ik);
            wfcpw->recip2real(tmpsi_in + this->max_npw, porter1, ik);
            std::complex<double> sup, sdown;
            if(this->veff->nc != 0)
            {
                const double* current_veff[4];
                for(int is=0;is<4;is++)
                {
                    current_veff[is] = &(this->veff[0](is, 0));
                }
                for (int ir = 0; ir < this->veff->nc; ir++)
                {
                    sup = porter[ir] * (current_veff[0][ir] + current_veff[3][ir])
                        + porter1[ir]
                                * (current_veff[1][ir]
                                - std::complex<double>(0.0, 1.0) * current_veff[2][ir]);
                    sdown = porter1[ir] * (current_veff[0][ir] - current_veff[3][ir])
                            + porter[ir]
                                * (current_veff[1][ir]
                                    + std::complex<double>(0.0, 1.0) * current_veff[2][ir]);
                    porter[ir] = sup;
                    porter1[ir] = sdown;
                }
            }
            // (3) fft back to G space.
            wfcpw->real2recip(porter, tmhpsi, this->ik, true);
            wfcpw->real2recip(porter1, tmhpsi + this->max_npw, this->ik, true);

            delete[] porter1;
        }
        tmhpsi += this->max_npw * this->npol;
        tmpsi_in += this->max_npw * this->npol;
    }
    delete[] porter;
    ModuleBase::timer::tick("Operator", "VeffPW");
    return;
}

} // namespace hamilt