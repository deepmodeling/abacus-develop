#include "spin_constrain.h"

template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::cal_mw_from_lambda()
{
    std::cout << "hello world" << std::endl;
    //this->pelec->charge->save_rho_before_sum_band();
    //this->phsol->solve(this->p_hamilt, this->psi[0], this->pelec, this->KS_SOLVER);
    //this->pelec->cal_energies(1);
    //this->pelec->f_en.deband = this->pelec->cal_delta_eband();
}

template class SpinConstrain<double, psi::DEVICE_CPU>;