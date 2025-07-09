#ifndef TD_FOLDING_H
#define TD_FOLDING_H

namespace module_rt{
// folding HR to hk, for hybrid gague
template<typename TR>
void folding_HR_td(const hamilt::HContainer<TR>& hR,
            std::complex<double>* hk,
            const ModuleBase::Vector3<double>& kvec_d_in,
            const int ncol,
            const int hk_type,
            const UnitCell& ucell,
            const ModuleBase::Vector3<double>& At);
}// namespace module_rt