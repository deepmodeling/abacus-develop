//=======================
// AUTHOR : jiyy
// DATE :   2024-01-10
//=======================

#ifndef AUXILIARY_FUNC_H
#define AUXILIARY_FUNC_H

#include "gaussian_abfs.h"
//#include "source_basis/module_pw/pw_basis_k.h"
#include "source_cell/klist.h"

#include <array>
#include <vector>

namespace Singular_Value
{
/**
     * @brief Calculating correction of Coulomb singularity
     * 
        carrier, // Phys. Rev. B, 75:205126, May 2007.
        massidda, // Phys. Rev. B 48, 5058. August 1993.
    */

using T_cal_fq_type = std::function<double(const ModuleBase::Vector3<double>& gk)>;
using T_cal_fq_type_no = std::function<double()>;

double cal_carrier(const UnitCell& ucell,
                  const std::vector<ModuleBase::Vector3<double>>& kvec_c,
                  const int& qdiv,
                  const double& qdense,
                  const int& niter,
                  const double& eps,
                  const int& a_rate);
double cal_massidda(const UnitCell& ucell,
                  const std::array<int, 3>& nmp,
                  const int& qdiv,
                  const double& start_lambda,
                  const int& niter,
                  const double& eps);

double solve_chi(const ModuleBase::Matrix3& G,
                 const std::vector<ModuleBase::Vector3<double>>& kvec_c,
                 const T_cal_fq_type& func_cal_fq,
                 const std::array<int, 3>& nq_arr,
                 const int& niter,
                 const double& eps,
                 const int& a_rate);
double solve_chi(const std::vector<ModuleBase::Vector3<double>>& kvec_c,
                 const T_cal_fq_type& func_cal_fq,
                 const double& fq_int);
double solve_chi(const int& nks, const T_cal_fq_type_no& func_cal_fq, const double& fq_int);
double sum_for_solve_chi(const std::vector<ModuleBase::Vector3<double>>& kvec_c,
                         const T_cal_fq_type& func_cal_fq,
                         const double& fq_int);
double Iter_Integral(const ModuleBase::Matrix3& G,
                     const T_cal_fq_type& func_cal_fq,
                     const std::array<int, 3>& nq_arr,
                     const int& niter,
                     const double& eps,
                     const int& a_rate);

// TODO: lower dimension please see PHYSICAL REVIEW B 87, 165122 (2013)

// qdiv=2 i.e. q^{-2} for 3D;
// qdiv=1 i.e. q^{-1} for 2D.
double fq_carrier(const double& tpiba,
                 const ModuleBase::Vector3<double>& qvec,
                 const int& qdiv,
                 std::vector<ModuleBase::Vector3<double>>& avec,
                 std::vector<ModuleBase::Vector3<double>>& bvec);
// gamma: chosen as the radius of sphere which has the same volume as the Brillouin zone.
double fq_massidda(const double& tpiba,
                 Gaussian_Abfs& gaussian_abfs,
                 const int& qdiv,
                 const double& lambda,
                 const int& lmax);
}; // namespace Singular_Value

#endif