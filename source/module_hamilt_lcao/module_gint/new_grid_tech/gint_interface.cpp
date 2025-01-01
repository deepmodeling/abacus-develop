#include "gint_interface.h"
#include "gint_vl.h"
#include "gint_vl_nspin4.h"
#include "gint_fvl.h"
#include "gint_rho.h"

namespace ModuleGint
{

void cal_gint_vl(
    const double* vr_eff,
    HContainer<double>* hR)
{
    Gint_vl gint_vl(vr_eff, hR);
    gint_vl.cal_gint();
}

void cal_gint_vl(
    std::vector<const double*> vr_eff,
    HContainer<std::complex<double>>* hR)
{
    Gint_vl_nspin4 gint_vl_nspin4(vr_eff, hR);
    gint_vl_nspin4.cal_gint();
}

void cal_gint_rho(
    const std::vector<HContainer<double>*>& DMR_vec,
    const int nspin,
    double **rho)
{
    Gint_rho gint_rho(DMR_vec, nspin, rho);
    gint_rho.cal_gint();
}

void cal_gint_fvl(
    const int nspin,
    const std::vector<const double*>& vr_eff,
    const std::vector<HContainer<double>*>& DMR_vec,
    const bool isforce,
    const bool isstress,
    ModuleBase::matrix* fvl,
    ModuleBase::matrix* svl)
{
    Gint_fvl gint_fvl(nspin, vr_eff, DMR_vec, isforce, isstress, fvl, svl);
    gint_fvl.cal_gint();
}

} // namespace ModuleGint