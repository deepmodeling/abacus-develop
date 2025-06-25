#pragma once

#include <memory>
#include <vector>
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "gint.h"
#include "gint_info.h"

namespace ModuleGint
{

class Gint_env_gamma : public Gint
{
    public:
    Gint_env_gamma(
        const double* psid,
        const Parallel_Orbitals* pv,
        const int nbands,
        const int nlocal,
        double* rho);

    void cal_env_band(const int iband);

    private:
    // output
    double* rho_;

    // intermediate variable
    std::vector<double> wfc_gint_;
};

}