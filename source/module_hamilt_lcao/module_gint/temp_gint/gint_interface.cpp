#include "gint_interface.h"
#include "module_base/timer.h"
#include "gint_vl.h"
#include "gint_vl_metagga.h"
#include "gint_vl_nspin4.h"
#include "gint_vl_metagga_nspin4.h"
#include "gint_fvl.h"
#include "gint_fvl_meta.h"
#include "gint_rho.h"
#include "gint_tau.h"

namespace ModuleGint
{

void cal_gint_vl(
    const double* vr_eff,
    HContainer<double>* hR)
{
    ModuleBase::timer::tick("Gint", "cal_gint_vl");
    Gint_vl gint_vl(vr_eff, hR);
    gint_vl.cal_gint();
    ModuleBase::timer::tick("Gint", "cal_gint_vl");
}

void cal_gint_vl(
    std::vector<const double*> vr_eff,
    HContainer<std::complex<double>>* hR)
{
    ModuleBase::timer::tick("Gint", "cal_gint_vl");
    Gint_vl_nspin4 gint_vl_nspin4(vr_eff, hR);
    gint_vl_nspin4.cal_gint();
    ModuleBase::timer::tick("Gint", "cal_gint_vl");
}

void cal_gint_vl_metagga(
    const double* vr_eff,
    const double* vfork,
    HContainer<double>* hR)
{
    ModuleBase::timer::tick("Gint", "cal_gint_vl_metagga");
    Gint_vl_metagga gint_vl_metagga(vr_eff, vfork, hR);
    gint_vl_metagga.cal_gint();
    ModuleBase::timer::tick("Gint", "cal_gint_vl_metagga");
}

void cal_gint_vl_metagga(
    std::vector<const double*> vr_eff,
    std::vector<const double*> vofk,
    HContainer<std::complex<double>>* hR)
{
    ModuleBase::timer::tick("Gint", "cal_gint_vl_metagga");
    Gint_vl_metagga_nspin4 gint_vl_metagga_nspin4(vr_eff, vofk, hR);
    gint_vl_metagga_nspin4.cal_gint();
    ModuleBase::timer::tick("Gint", "cal_gint_vl_metagga");
}

void cal_gint_rho(
    const std::vector<HContainer<double>*>& dm_vec,
    const int nspin,
    double **rho)
{
    ModuleBase::timer::tick("Gint", "cal_gint_rho");
    Gint_rho gint_rho(dm_vec, nspin, rho);
    gint_rho.cal_gint();
    ModuleBase::timer::tick("Gint", "cal_gint_rho");
}

void cal_gint_tau(        
    const std::vector<HContainer<double>*>& dm_vec,
    const int nspin,
    double** tau)
{
    ModuleBase::timer::tick("Gint", "cal_gint_tau");
    Gint_tau gint_tau(dm_vec, nspin, tau);
    gint_tau.cal_gint();
    ModuleBase::timer::tick("Gint", "cal_gint_tau");
}

void cal_gint_fvl(
    const int nspin,
    const std::vector<const double*>& vr_eff,
    const std::vector<HContainer<double>*>& dm_vec,
    const bool isforce,
    const bool isstress,
    ModuleBase::matrix* fvl,
    ModuleBase::matrix* svl)
{
    ModuleBase::timer::tick("Gint", "cal_gint_fvl");
    Gint_fvl gint_fvl(nspin, vr_eff, dm_vec, isforce, isstress, fvl, svl);
    gint_fvl.cal_gint();
    ModuleBase::timer::tick("Gint", "cal_gint_fvl");
}

void cal_gint_fvl_meta(
    const int nspin,
    const std::vector<const double*>& vr_eff,
    const std::vector<const double*>& vofk,
    const std::vector<HContainer<double>*>& dm_vec,
    const bool isforce,
    const bool isstress,
    ModuleBase::matrix* fvl,
    ModuleBase::matrix* svl)
{
    ModuleBase::timer::tick("Gint", "cal_gint_fvl_meta");
    Gint_fvl_meta gint_fvl_meta(nspin, vr_eff, vofk, dm_vec, isforce, isstress, fvl, svl);
    gint_fvl_meta.cal_gint();
    ModuleBase::timer::tick("Gint", "cal_gint_fvl_meta");
}

} // namespace ModuleGint