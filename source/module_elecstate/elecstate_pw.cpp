#include "elecstate_pw.h"

namespace elecstate
{

void ElecStatePW::init(Charge* chg_in)
{
    this->pchg = chg_in;
}

const hamilt::MatrixBlock<double> ElecStatePW::getRho() const
{
    hamilt::MatrixBlock<double> temp{&(this->pchg->rho[0][0]), 1, 1}; // this->chr->get_nspin(), this->chr->get_nrxx()};
    return temp;
}

void ElecStatePW::updateRhoK(const psi::Psi<std::complex<double>>& psi)
{
    this->pchg->sum_band_k(/*psi*/);
    return;
}

void ElecStatePW::getNewRho()
{
    return;
}

} // namespace elecstate