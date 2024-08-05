#include "gint_inout.h"

Gint_inout::Gint_inout(double** rho_in, 
                       Gint_Tools::job_type job_in, 
                    bool if_symm_in )
{
    rho = rho_in;
    job = job_in;
    if_symm = if_symm_in;
}

Gint_inout::Gint_inout(const int ispin_in,
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

Gint_inout::Gint_inout(const int ispin_in,
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

Gint_inout::Gint_inout(const double* vl_in, int ispin_in, Gint_Tools::job_type job_in)
{
    vl = vl_in;
    ispin = ispin_in;
    job = job_in;
}

Gint_inout::Gint_inout(const double* vl_in, 
            const double* vofk_in, 
            int ispin_in, 
            Gint_Tools::job_type job_in)
{
    vl = vl_in;
    vofk = vofk_in;
    ispin = ispin_in;
    job = job_in;
}

Gint_inout::Gint_inout(double*** DM_in,
              double** rho_in, 
              Gint_Tools::job_type job_in, 
              bool if_symm_in)
{
    DM = DM_in;
    rho = rho_in;
    job = job_in;
    if_symm = if_symm_in;
}

Gint_inout::Gint_inout(const double* vl_in, 
                        Gint_Tools::job_type job_in)
{
    vl = vl_in;
    job = job_in;
}

Gint_inout::Gint_inout(const double* vl_in, 
                        const double* vofk_in, 
                        Gint_Tools::job_type job_in)
{
    vl = vl_in;
    vofk = vofk_in;
    job = job_in;
}