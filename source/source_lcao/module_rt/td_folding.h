#ifndef TD_FOLDING_H
#define TD_FOLDING_H
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "source_base/abfs-vector3_order.h"

namespace module_rt{
// folding HR to hk, for hybrid gauge
template<typename TR>
void folding_HR_td(const UnitCell& ucell,
            const hamilt::HContainer<TR>& hR,
            std::complex<double>* hk,
            const ModuleBase::Vector3<double>& kvec_d_in,
            const ModuleBase::Vector3<double>& At,
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
            const int ix,
            const int ncol,
            const int hk_type);
}// namespace module_rt

#endif