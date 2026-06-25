#ifndef TD_FOLDING_H
#define TD_FOLDING_H
#include "source_lcao/module_hcontainer/hcontainer.h"

namespace module_rt{
// folding HR to hk, for hybrid gauge
template<typename TR>
void folding_HR_td(const hamilt::HContainer<TR>& hR,
                std::complex<double>* hk,
                const ModuleBase::Vector3<double>& kvec_d_in,
                const ModuleBase::Vector3<double>& At,
                const std::map<ModuleBase::Vector3<int>, std::complex<double>>& phase_hybrid,
                const int ncol,
                const int hk_type);
template<typename TR>
void folding_partial_HR(const UnitCell& ucell,
                const hamilt::HContainer<TR>& hR,
                std::complex<double>* hk,
                const ModuleBase::Vector3<double>& kvec_d_in,
                const int ix,
                const int ncol,
                const int hk_type);
template<typename TR>
void folding_partial_HR_td(const UnitCell& ucell,
            const hamilt::HContainer<TR>& hR,
            std::complex<double>* hk,
            const ModuleBase::Vector3<double>& kvec_d_in,
            const ModuleBase::Vector3<double>& cart_At,
            const std::map<ModuleBase::Vector3<int>, std::complex<double>>& phase_hybrid,
            const int ix,
            const int ncol,
            const int hk_type);
void folding_partial_dot(const hamilt::HContainer<double>& dR,
            std::complex<double>* dk,
            const ModuleBase::Vector3<double>& kvec_d_in,
            const int ncol,
            const int hk_type,
            const UnitCell* ucell,
            const std::map<ModuleBase::Vector3<int>, std::complex<double>>& phase_hybrid,
            const ModuleBase::Vector3<double>& At,
            const ModuleBase::Vector3<double>& Et);
}// namespace module_rt

#endif