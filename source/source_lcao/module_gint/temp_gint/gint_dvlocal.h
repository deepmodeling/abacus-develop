#pragma once
#include <memory>
#include <vector>
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "source_lcao/hamilt_lcaodft/LCAO_HS_arrays.hpp"
#include "source_base/abfs-vector3_order.h"
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
        const int npol)
        : vr_eff_(vr_eff), nspin_(nspin), npol_(npol), dr3_(gint_info_->get_mgrid_volume())
        {
            assert(nspin_ == 2); //   currently only npin == 2 is supported
        }
    
    void cal_dvlocal();

    void cal_dvlocal_R_sparseMatrix(
        const int nspin,
        const int cspin,
        const int nlocal,
        const double sparse_thr, 
        const Parallel_Orbitals& pv,
        const UnitCell& ucell,
        const Grid_Driver& gdriver,
        LCAO_HS_Arrays& hs_arrays);
    
    private:
    void init_hr_gint_();

    void cal_hr_gint_();

    void distribute_pvdpR_sparseMatrix(
        const int cspin,
        const int dim,
        const int nlocal,
        const double sparse_threshold,
        const std::map<Abfs::Vector3_Order<int>,
                       std::map<size_t, std::map<size_t, double>>>&
            pvdpR_sparseMatrix,
        const Parallel_Orbitals& pv,
        LCAO_HS_Arrays& HS_Arrays);

    // input
    const double* vr_eff_;
    int nspin_;
    int npol_;

    // intermediate variables
    double dr3_;
    HContainer<double> pvdpRx;
    HContainer<double> pvdpRy;
    HContainer<double> pvdpRz;
};

}