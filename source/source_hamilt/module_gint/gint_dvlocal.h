#ifndef ABACUS_SOURCE_HAMILT_MODULE_GINT_GINT_DVLOCAL_H
#define ABACUS_SOURCE_HAMILT_MODULE_GINT_GINT_DVLOCAL_H

#include <memory>
#include <vector>
#include "source_hamilt/module_hcontainer/hcontainer.h"
#include "gint.h"
#include "gint_info.h"

namespace ModuleGint
{

class Gint_dvlocal : public Gint
{
    public:
    Gint_dvlocal(
        const double* vr_eff,
        const int nspin,
        const int npol,
        const bool full_triangle = false)
        : vr_eff_(vr_eff), nspin_(nspin), npol_(npol), full_triangle_(full_triangle),
          dr3_(gint_info_->get_mgrid_volume())
        {
            assert(nspin_ == 1 || nspin_ == 2); //   currently only nspin == 1 or 2 is supported
        }
    
    void cal_dvlocal();

    HContainer<double>* get_pvdpRx() { return &pvdpRx; }
    HContainer<double>* get_pvdpRy() { return &pvdpRy; }
    HContainer<double>* get_pvdpRz() { return &pvdpRz; }
    
    private:
    void init_hr_gint_();

    void cal_hr_gint_();

    // input
    const double* vr_eff_ = nullptr;
    int nspin_;
    int npol_;
    // if true, fill both triangles of pvdpR (all directed atom pairs); default upper-only
    bool full_triangle_ = false;

    // intermediate variables
    double dr3_;
    HContainer<double> pvdpRx;
    HContainer<double> pvdpRy;
    HContainer<double> pvdpRz;
};

}

#endif // ABACUS_SOURCE_HAMILT_MODULE_GINT_GINT_DVLOCAL_H
