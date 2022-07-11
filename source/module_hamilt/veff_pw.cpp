#include "veff_pw.h"

#include "module_base/timer.h"
#include "src_pw/global.h"
#include "module_base/tool_quit.h"

namespace hamilt
{

VeffPW::VeffPW(
    int max_npw_in,
    int npol_in,
    const int* ngk_in,
    const int* isk_in,
    const ModuleBase::matrix* veff_in,
    ModulePW::PW_Basis_K* wfcpw_in
)
{
    this->max_npw = max_npw_in;
    this->npol = npol_in;
    this->ngk = ngk_in;
    this->isk = isk_in;
    this->veff = veff_in;
    this->wfcpw = wfcpw_in;
    if(this->max_npw == 0 || this->npol == 0 || this->ngk == nullptr
    || this->isk == nullptr || this->veff == nullptr || this->wfcpw == nullptr)
    {
        ModuleBase::WARNING_QUIT("VeffPW", "Constuctor of Operator::VeffPW is failed, please check your code!");
    }
}

void VeffPW::act(const std::complex<double> *psi_in, std::complex<double> *hpsi, const size_t size) const
{
    ModuleBase::timer::tick("Operator", "VeffPW");
    int m = int(size / this->max_npw / this->npol);
    if (int(size - m * this->max_npw * this->npol) != 0)
    {
        m++;
    }
    const int npw = this->ngk[this->ik];

    const int current_spin = this->isk[this->ik];

    std::complex<double> *tmhpsi;
    const std::complex<double> *tmpsi_in;

    std::complex<double> *porter = new std::complex<double>[wfcpw->nmaxgr];
    tmhpsi = hpsi;
    tmpsi_in = psi_in;
    for (int ib = 0; ib < m; ++ib)
    {
        if (this->npol == 1)
        {
            const double* current_veff = &(this->veff[0](current_spin, 0));
            wfcpw->recip2real(tmpsi_in, porter, ik);
            for (int ir = 0; ir < this->veff->nc; ++ir)
            {
                porter[ir] *= current_veff[ir];
            }
            wfcpw->real2recip(porter, tmhpsi, ik, true);
        }
        else
        {
            const double* current_veff[4];
            for(int is=0;is<4;is++)
            {
                current_veff[is] = &(this->veff[0](is, 0));
            }
            std::complex<double> *porter1 = new std::complex<double>[wfcpw->nmaxgr];
            // fft to real space and doing things.
            wfcpw->recip2real(tmpsi_in, porter, ik);
            wfcpw->recip2real(tmpsi_in + this->max_npw, porter1, ik);
            std::complex<double> sup, sdown;
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