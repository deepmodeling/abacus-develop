#include "hamilt_pw.h"

namespace hamilt
{

void HamiltPW::ch_mock()
{
    return;
}

void HamiltPW::hk_mock()
{
    return;
}

void HamiltPW::hpsi_mock(const psi::Psi<std::complex<double>>& psi, psi::Psi<std::complex<double>>& hpsi) const
{
    for (size_t iband = 0; iband < hpsi.get_nbands(); ++iband)
    {
        for (size_t ibasis = 0; ibasis < hpsi.get_nbasis(); ++ibasis)
        {
            hpsi(iband, ibasis) = psi(iband, ibasis);
        }
    }
}

} // namespace hamilt