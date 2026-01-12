#pragma once

#include <memory>
#include <vector>
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "gint.h"
#include "gint_info.h"

namespace ModuleGint
{

class Gint_env_k : public Gint
{
    public:
    Gint_env_k(
        const std::complex<double>* psid,
        const Parallel_Orbitals* pv,
        const std::vector<Vec3d>& kvec_c,
        const std::vector<Vec3d>& kvec_d,
        const int nbands,
        const int nlocal,
        const int ik,
        const int nspin,
        const int npol,
        double* rho);

    void cal_env_band(const int iband);

    private:
    // input
    const std::vector<Vec3d>& kvec_c_;
    const std::vector<Vec3d>& kvec_d_;
    int ik_;
    int nspin_;
    int npol_;

    // output
    double* rho_;

    // intermediate variable
    std::vector<std::complex<double>> wfc_gint_;
};

}