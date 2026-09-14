#pragma once

#include "gint.h"
#include "gint_info.h"

#include <memory>
#include <vector>

namespace ModuleGint
{

class Gint_env_k : public Gint
{
  public:
    Gint_env_k(const std::complex<double>* psid,
               const Parallel_Orbitals* pv,
               const std::vector<Vec3d>& kvec_d,
               const int nbands,
               const int nlocal,
               const int ik,
               const int npol,
               std::complex<double>* wfc);

    void cal_env_band(const int iband);

  private:
    // input
    const std::vector<Vec3d>& kvec_d_;
    int ik_;
    int npol_;

    // output
    std::complex<double>* wfc_ = nullptr;

    // intermediate variable
    std::vector<std::complex<double>> wfc_gint_;
};

} // namespace ModuleGint
