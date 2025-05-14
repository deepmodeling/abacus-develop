//=======================
// AUTHOR : jiyy
// DATE :   2024-01-10
//=======================

#ifndef AUXILIARY_FUNC_H
#define AUXILIARY_FUNC_H

#include "gaussian_abfs.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_cell/klist.h"

#include <array>
#include <vector>

class Singular_Value
{
  public:
    enum class Fq_type
    {
        Type_0, // Phys. Rev. B, 75:205126, May 2007.
        Type_1, // Phys. Rev. B 48, 5058. August 1993.
    };

  private:
    using T_cal_fq_type = std::function<double(const ModuleBase::Vector3<double>& gk)>;
    using T_cal_fq_type_no = std::function<double()>;

  public:
    static double cal_type_0(const UnitCell& ucell,
                             const std::vector<ModuleBase::Vector3<double>>& kvec_c,
                             const int& qdiv,
                             const double& qdense,
                             const int& niter,
                             const double& eps,
                             const int& a_rate);
    static double cal_type_1(const UnitCell& ucell,
                             const std::array<int, 3>& nmp,
                             const int& qdiv,
                             const double& start_lambda,
                             const int& niter,
                             const double& eps);

  private:
    static double solve_chi(const ModuleBase::Matrix3 &G,
                            const std::vector<ModuleBase::Vector3<double>>& kvec_c,
                            const T_cal_fq_type& func_cal_fq,
                            const std::array<int, 3>& nq_arr,
                            const int& niter,
                            const double& eps,
                            const int& a_rate);
    static double solve_chi(const std::vector<ModuleBase::Vector3<double>>& kvec_c,
                            const T_cal_fq_type& func_cal_fq,
                            const double& fq_int);
    static double solve_chi(const int& nks, const T_cal_fq_type_no& func_cal_fq, const double& fq_int);
    static double sum_for_solve_chi(const std::vector<ModuleBase::Vector3<double>>& kvec_c,
                                    const T_cal_fq_type& func_cal_fq,
                                    const double& fq_int);
    static double Iter_Integral(const ModuleBase::Matrix3 &G,
                                const T_cal_fq_type& func_cal_fq,
                                const std::array<int, 3>& nq_arr,
                                const int& niter,
                                const double& eps,
                                const int& a_rate);

    // TODO: lower dimension please see PHYSICAL REVIEW B 87, 165122 (2013)

    // qdiv=2 i.e. q^{-2} for 3D;
    // qdiv=1 i.e. q^{-1} for 2D.
    static double fq_type_0(const double& tpiba,
                            const ModuleBase::Vector3<double>& qvec,
                            const int& qdiv,
                            std::vector<ModuleBase::Vector3<double>>& avec,
                            std::vector<ModuleBase::Vector3<double>>& bvec);
    // gamma: chosen as the radius of sphere which has the same volume as the Brillouin zone.
    static double fq_type_1(const double& tpiba, Gaussian_Abfs& gaussian_abfs, const int& qdiv, const double& lambda, const int& lmax);
};

#endif