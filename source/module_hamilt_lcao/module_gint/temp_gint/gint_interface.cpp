#include "gint_interface.h"
#include "source_base/timer.h"
#include "module_parameter/parameter.h"
#include "gint_vl.h"
#include "gint_vl_metagga.h"
#include "gint_vl_nspin4.h"
#include "gint_vl_metagga_nspin4.h"
#include "gint_fvl.h"
#include "gint_fvl_meta.h"
#include "gint_rho.h"
#include "gint_tau.h"

#ifdef __CUDA
#include "gint_vl_gpu.h"
#include "gint_rho_gpu.h"
#include "gint_fvl_gpu.h"
#include "gint_vl_nspin4_gpu.h"
#include "gint_vl_metagga_gpu.h"
#include "gint_vl_metagga_nspin4_gpu.h"
#include "gint_tau_gpu.h"
#include "gint_fvl_meta_gpu.h"
#endif

namespace ModuleGint
{

void cal_gint_vl(
    const double* vr_eff,
    HContainer<double>* hR)
{
#ifdef __CUDA
    if(PARAM.inp.device == "gpu")
    {
        Gint_vl_gpu gint_vl(vr_eff, hR);
        gint_vl.cal_gint();
    } else
#endif
    {
        Gint_vl gint_vl(vr_eff, hR);
        gint_vl.cal_gint();
    }
}

// nspin == 4 case
void cal_gint_vl(
    std::vector<const double*> vr_eff,
    HContainer<std::complex<double>>* hR)
{
    #ifdef __CUDA
    if(PARAM.inp.device == "gpu")
    {
        Gint_vl_nspin4_gpu gint_vl_nspin4(vr_eff, hR);
        gint_vl_nspin4.cal_gint();
    } else
    #endif
    {
        Gint_vl_nspin4 gint_vl_nspin4(vr_eff, hR);
        gint_vl_nspin4.cal_gint();
    }
}

void cal_gint_vl_metagga(
    const double* vr_eff,
    const double* vfork,
    HContainer<double>* hR)
{
#ifdef __CUDA
    if(PARAM.inp.device == "gpu")
    {
        Gint_vl_metagga_gpu gint_vl_metagga(vr_eff, vfork, hR);
        gint_vl_metagga.cal_gint();
    } else
#endif
    {
        Gint_vl_metagga gint_vl_metagga(vr_eff, vfork, hR);
        gint_vl_metagga.cal_gint();
    }
}

// nspin == 4 case
void cal_gint_vl_metagga(
    std::vector<const double*> vr_eff,
    std::vector<const double*> vofk,
    HContainer<std::complex<double>>* hR)
{
#ifdef __CUDA
    if(PARAM.inp.device == "gpu")
    {
        Gint_vl_metagga_nspin4_gpu gint_vl_metagga_nspin4(vr_eff, vofk, hR);
        gint_vl_metagga_nspin4.cal_gint();
    } else
#endif
    {
        Gint_vl_metagga_nspin4 gint_vl_metagga_nspin4(vr_eff, vofk, hR);
        gint_vl_metagga_nspin4.cal_gint();
    }
}

void cal_gint_rho(
    const std::vector<HContainer<double>*>& dm_vec,
    const int nspin,
    double **rho,
    bool is_dm_symm)
{
    #ifdef __CUDA
    if(PARAM.inp.device == "gpu")
    {
        Gint_rho_gpu gint_rho(dm_vec, nspin, rho, is_dm_symm);
        gint_rho.cal_gint();
    } else
    #endif
    {
        Gint_rho gint_rho(dm_vec, nspin, rho, is_dm_symm);
        gint_rho.cal_gint();
    }
}

void cal_gint_tau(        
    const std::vector<HContainer<double>*>& dm_vec,
    const int nspin,
    double** tau)
{
    #ifdef __CUDA
    if(PARAM.inp.device == "gpu")
    {
        Gint_tau_gpu gint_tau(dm_vec, nspin, tau);
        gint_tau.cal_gint();
    } else
    #endif
    {
        Gint_tau gint_tau(dm_vec, nspin, tau);
        gint_tau.cal_gint();
    }
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
#ifdef __CUDA
    if(PARAM.inp.device == "gpu")
    {
        Gint_fvl_gpu gint_fvl_gpu(nspin, vr_eff, dm_vec, isforce, isstress, fvl, svl);
        gint_fvl_gpu.cal_gint();
    } else
#endif
    {
        Gint_fvl gint_fvl(nspin, vr_eff, dm_vec, isforce, isstress, fvl, svl);
        gint_fvl.cal_gint();
    }
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
#ifdef __CUDA
    if(PARAM.inp.device == "gpu")
    {
        Gint_fvl_meta_gpu gint_fvl_meta(nspin, vr_eff, vofk, dm_vec, isforce, isstress, fvl, svl);
        gint_fvl_meta.cal_gint();
    } else
#endif
    {
        Gint_fvl_meta gint_fvl_meta(nspin, vr_eff, vofk, dm_vec, isforce, isstress, fvl, svl);
        gint_fvl_meta.cal_gint();
    }
}

void cal_dvlocal_R_sparseMatrix(
    const int nspin,
    const int npol,
    const int current_spin,
    const int nlocal,
    const double sparse_thr,
    const double* vr_eff,
    const Parallel_Orbitals& pv,
    const UnitCell& ucell,
    const Grid_Driver& gdriver,
    LCAO_HS_Arrays& hs_arrays)
{
    Gint_dvlocal gint_dvlocal(vr_eff, nspin, npol);
    gint_dvlocal.cal_dvlocal();
    gint_dvlocal.cal_dvlocal_R_sparseMatrix(
        nspin, current_spin, nlocal, sparse_thr,
        pv, ucell, gdriver, hs_arrays);
}

} // namespace ModuleGint