#include "elecstatepw.h"

namespace ModuleElecS
{

void ElecStatePW::init(
    Charge* chg_in
    )
{
    this->pchg = chg_in;
}

const ModuleHamilt::MatrixBlock<double> ElecStatePW::getRho()const
{
    ModuleHamilt::MatrixBlock<double> temp{&(this->pchg->rho[0][0]), 1,1};//this->chr->get_nspin(), this->chr->get_nrxx()};
    return temp;
}

void ElecStatePW::updateRhoK(const ModulePsi::Psi<std::complex<double>> &psi)
{
    this->pchg->sum_band_k(/*psi*/);
    return;
}

void ElecStatePW::getNewRho()
{
    return;
}

}