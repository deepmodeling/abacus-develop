#ifndef GINT_INOUT_HPP
#define GINT_INOUT_HPP
#include "gint_tools.h"
#include "module_elecstate/module_charge/charge.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_base/array_pool.h"

#include <cstdlib>
// the class is used to pass input/output variables
// into the unified interface gint
// not sure if this is the best practice though ..
class Gint_inout
{
  public:
    // input
    double*** DM;
    const double* vl;
    const double* vofk;
    bool isforce;
    bool isstress;
    int ispin;
    bool if_symm = false; // if true, use dsymv in gint_kernel_rho; if false, use dgemv.

    // output
    double** rho;
    ModuleBase::matrix* fvl_dphi;
    ModuleBase::matrix* svl_dphi;
    Gint_Tools::job_type job;

    // electron density and kin_r, multi-k
    Gint_inout(double** rho_in, Gint_Tools::job_type job_in, bool if_symm_in = true)
    {
        rho = rho_in;
        job = job_in;
        if_symm = if_symm_in;
    }

    // force
    Gint_inout(const int ispin_in,
               const double* vl_in,
               bool isforce_in,
               bool isstress_in,
               ModuleBase::matrix* fvl_dphi_in,
               ModuleBase::matrix* svl_dphi_in,
               Gint_Tools::job_type job_in)
    {
        vl = vl_in;
        isforce = isforce_in;
        isstress = isstress_in;
        fvl_dphi = fvl_dphi_in;
        svl_dphi = svl_dphi_in;
        job = job_in;
        ispin = ispin_in;
    }

    // force (mGGA)
    Gint_inout(const int ispin_in,
               const double* vl_in,
               const double* vofk_in,
               const bool isforce_in,
               const bool isstress_in,
               ModuleBase::matrix* fvl_dphi_in,
               ModuleBase::matrix* svl_dphi_in,
               Gint_Tools::job_type job_in)
    {
        vl = vl_in;
        vofk = vofk_in;
        isforce = isforce_in;
        isstress = isstress_in;
        fvl_dphi = fvl_dphi_in;
        svl_dphi = svl_dphi_in;
        job = job_in;
        ispin = ispin_in;
    }

    // vlocal, multi-k
    Gint_inout(const double* vl_in, int ispin_in, Gint_Tools::job_type job_in)
    {
        vl = vl_in;
        ispin = ispin_in;
        job = job_in;
    }

    // mGGA vlocal, multi-k
    Gint_inout(const double* vl_in, const double* vofk_in, int ispin_in, Gint_Tools::job_type job_in)
    {
        vl = vl_in;
        vofk = vofk_in;
        ispin = ispin_in;
        job = job_in;
    }

    // electron density and kin_r, gamma point
    Gint_inout(double*** DM_in, double** rho_in, Gint_Tools::job_type job_in, bool if_symm_in = true)
    {
        DM = DM_in;
        rho = rho_in;
        job = job_in;
        if_symm = if_symm_in;
    }

    // vlocal, gamma point
    Gint_inout(const double* vl_in, Gint_Tools::job_type job_in)
    {
        vl = vl_in;
        job = job_in;
    }

    // mGGA vlocal, gamma point
    Gint_inout(const double* vl_in, const double* vofk_in, Gint_Tools::job_type job_in)
    {
        vl = vl_in;
        vofk = vofk_in;
        job = job_in;
    }
};
#endif