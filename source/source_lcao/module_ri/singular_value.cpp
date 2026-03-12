//=======================
// AUTHOR : jiyy
// DATE :   2024-01-01
//=======================

#ifndef AUXILIARY_FUNC_CPP
#define AUXILIARY_FUNC_CPP

#include "singular_value.h"

#include "conv_coulomb_pot_k.h"
#include "source_base/global_variable.h"
#include "source_base/math_ylmreal.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
//#include "source_pw/hamilt_pwdft/global.h"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace Singular_Value
{
// for analytic integral of fq
double sum_for_solve_chi(const std::vector<ModuleBase::Vector3<double>>& kvec_c,
                         const T_cal_fq_type& func_cal_fq,
                         const double& fq_int)
{
    const int nks = kvec_c.size();

    // cal fq sum except q=0
    double fq_sum = 0;
    for (size_t ik = 0; ik != nks; ++ik)
        fq_sum += kvec_c[ik].norm() ? func_cal_fq(kvec_c[ik]) : 0;

    double chi = fq_int * nks - fq_sum;

    return chi;
}

// for numerical integral of fq
double solve_chi(const ModuleBase::Matrix3& G,
                 const std::vector<ModuleBase::Vector3<double>>& kvec_c,
                 const T_cal_fq_type& func_cal_fq,
                 const std::array<int, 3>& nq_arr,
                 const int& niter,
                 const double& eps,
                 const int& a_rate)
{
    // cal fq integral
    double fq_int = Iter_Integral(G, func_cal_fq, nq_arr, niter, eps, a_rate);

    return sum_for_solve_chi(kvec_c, func_cal_fq, fq_int);
}

// for analytic integral of fq
double solve_chi(const std::vector<ModuleBase::Vector3<double>>& kvec_c,
                 const T_cal_fq_type& func_cal_fq,
                 const double& fq_int)
{
    return sum_for_solve_chi(kvec_c, func_cal_fq, fq_int);
}

// for analytic integral of fq with gaussian sum
double solve_chi(const int& nks, const T_cal_fq_type_no& func_cal_fq, const double& fq_int)
{
    double chi = fq_int * nks - func_cal_fq();

    return chi;
}

double fq_carrier(const double& tpiba,
                 const ModuleBase::Vector3<double>& qvec,
                 const int& qdiv,
                 std::vector<ModuleBase::Vector3<double>>& avec,
                 std::vector<ModuleBase::Vector3<double>>& bvec)
{
    assert(qvec.norm2());

    std::vector<ModuleBase::Vector3<double>> baq(3);
    std::vector<ModuleBase::Vector3<double>> baq_2(3);
    const int qexpo = -abs(qdiv);
    const double prefactor = std::pow(ModuleBase::TWO_PI, -qexpo);

    for (size_t i = 0; i != 3; ++i)
    {
        baq[i] = tpiba * bvec[i] * std::sin(avec[i] * qvec * ModuleBase::TWO_PI);
        baq_2[i] = tpiba * bvec[i] * std::sin(avec[i] * qvec * ModuleBase::PI);
    }

    double sum_baq = 0;
    double sum_baq_2 = 0;
    for (size_t i = 1; i != 4; ++i)
    {
        size_t j = i % 3 + 1;
        size_t new_i = i - 1;
        size_t new_j = j - 1;
        sum_baq_2 += baq_2[new_i] * baq_2[new_i];
        sum_baq += baq[new_i] * baq[new_j];
    }
    double fq = prefactor * std::pow(4 * sum_baq_2 + 2 * sum_baq, 0.5 * qexpo);

    return fq;
}

double cal_carrier(const UnitCell& ucell,
                  const std::vector<ModuleBase::Vector3<double>>& kvec_c,
                  const int& qdiv,
                  const double& qdense,
                  const int& niter,
                  const double& eps,
                  const int& a_rate)
{
    ModuleBase::TITLE("Singular_Value", "cal_carrier");
    ModuleBase::timer::tick("Singular_Value", "cal_carrier");

    std::vector<ModuleBase::Vector3<double>> avec = {ucell.a1, ucell.a2, ucell.a3};
    std::vector<ModuleBase::Vector3<double>> bvec;
    bvec.resize(3);
    bvec[0].x = ucell.G.e11;
    bvec[0].y = ucell.G.e12;
    bvec[0].z = ucell.G.e13;

    bvec[1].x = ucell.G.e21;
    bvec[1].y = ucell.G.e22;
    bvec[1].z = ucell.G.e23;

    bvec[2].x = ucell.G.e31;
    bvec[2].y = ucell.G.e32;
    bvec[2].z = ucell.G.e33;

    std::array<int, 3> nq_arr;
    const double qdense_tpiba = qdense * ucell.tpiba;
    std::transform(bvec.begin(),
                   bvec.end(),
                   nq_arr.begin(),
                   [&qdense_tpiba, &a_rate](ModuleBase::Vector3<double>& vec) -> int {
                       int index = static_cast<int>(vec.norm() * qdense_tpiba);
                       return index ? index - index % a_rate : a_rate;
                   });
    const T_cal_fq_type func_cal_fq_carrier
        = std::bind(&fq_carrier, ucell.tpiba, std::placeholders::_1, qdiv, avec, bvec);

    double val = solve_chi(ucell.G, kvec_c, func_cal_fq_carrier, nq_arr, niter, eps, a_rate);
    ModuleBase::timer::tick("Singular_Value", "cal_carrier");
    return val;
}

double fq_massidda(const double& tpiba,
                 Gaussian_Abfs& gaussian_abfs,
                 const int& qdiv,
                 const double& lambda,
                 const int& lmax)
{
    const size_t ik = 0;
    const double qexpo = -abs(qdiv);
    const bool exclude_Gamma = true;
    const ModuleBase::Vector3<double> tau(0, 0, 0);

    auto lattice_sum = gaussian_abfs.get_lattice_sum(tpiba, ik, qexpo, lambda, exclude_Gamma, lmax, tau);
    assert(lattice_sum[0].imag() < 1e-10);
    double fq = lattice_sum[0].real() * std::sqrt(ModuleBase::FOUR_PI);

    return fq;
}

double cal_massidda(const UnitCell& ucell,
                  const std::array<int, 3>& nmp,
                  const int& qdiv,
                  const double& start_lambda,
                  const int& niter,
                  const double& eps)
{
    ModuleBase::TITLE("Singular_Value", "cal_massidda");
    ModuleBase::timer::tick("Singular_Value", "cal_massidda");

    ModuleBase::Matrix3 bvec;
    bvec.e11 = ucell.G.e11 / nmp[0];
    bvec.e12 = ucell.G.e12 / nmp[0];
    bvec.e13 = ucell.G.e13 / nmp[0];

    bvec.e21 = ucell.G.e21 / nmp[1];
    bvec.e22 = ucell.G.e22 / nmp[1];
    bvec.e23 = ucell.G.e23 / nmp[1];

    bvec.e31 = ucell.G.e31 / nmp[2];
    bvec.e32 = ucell.G.e32 / nmp[2];
    bvec.e33 = ucell.G.e33 / nmp[2];

    const int nks = nmp[0] * nmp[1] * nmp[2];
    const std::vector<ModuleBase::Vector3<double>> qvec(1, ModuleBase::Vector3<double>{0, 0, 0});
    const int lmax = 0;

    auto cal_chi = [&qdiv, &bvec, &nks, &qvec, &lmax, &ucell](const double& lambda) {
        Gaussian_Abfs gaussian_abfs;
        const double exponent = 1 / lambda;
        gaussian_abfs.init(ucell, lmax, qvec, bvec, exponent);
        const T_cal_fq_type_no func_cal_fq_massidda
            = std::bind(&fq_massidda, ucell.tpiba, gaussian_abfs, qdiv, lambda, lmax);
        double prefactor
            = ModuleBase::TWO_PI * std::pow(lambda, -1.0 / qdiv) * ucell.omega / std::pow(ModuleBase::TWO_PI, 3);
        double fq_int;
        if (qdiv == 2)
            fq_int = prefactor * std::sqrt(ModuleBase::PI);
        else if (qdiv == 1)
            fq_int = prefactor;
        else
            ModuleBase::WARNING_QUIT("Singular_Value::cal_massidda", "Type 1 fq only supports qdiv=1 or qdiv=2!");
        return solve_chi(nks, func_cal_fq_massidda, fq_int);
    };

    int tot_iter = 0;
    double val_extra_old = 0.5 * std::numeric_limits<double>::max();
    double lammda_old = start_lambda;
    double val_old = cal_chi(lammda_old);
    double val_extra;
    for (size_t iter = 0; iter != niter; ++iter)
    {
        double lammda_new = lammda_old * 0.5;
        double val_new = cal_chi(lammda_new);
        double dval = (val_new - val_old) / (lammda_new - lammda_old);
        val_extra = val_new + dval * (0.0 - lammda_new);
        if (std::abs(val_extra - val_extra_old) < eps)
            break;
        lammda_old = lammda_new;
        val_old = val_new;
        val_extra_old = val_extra;
        tot_iter += 1;
    }

    if (tot_iter == niter)
        ModuleBase::WARNING_QUIT("Singular_Value::cal_massidda", "not converged!");

    ModuleBase::timer::tick("Singular_Value", "cal_massidda");
    return val_extra;
}

double Iter_Integral(const ModuleBase::Matrix3& G,
                     const T_cal_fq_type& func_cal_fq,
                     const std::array<int, 3>& nq_arr,
                     const int& niter,
                     const double& eps,
                     const int& a_rate)
{
    bool any_negative = std::any_of(nq_arr.begin(), nq_arr.end(), [](int i) { return i < 0; });
    bool any_nthree = std::any_of(nq_arr.begin(), nq_arr.end(), [&a_rate](int i) { return i % a_rate != 0; });
    if (any_negative || any_nthree)
        ModuleBase::WARNING_QUIT("Singular_Value::Iter_Integral",
                                 "The elements of `nq_arr` should be "
                                 "non-negative and multiples of a_rate!");
    bool all_zero = std::all_of(nq_arr.begin(), nq_arr.end(), [](int i) { return i == 0; });
    if (all_zero)
        ModuleBase::WARNING_QUIT("Singular_Value::Iter_Integral",
                                 "At least one element of `nq_arr` should be non-zero!");

    const int nqs = std::accumulate(nq_arr.begin(), nq_arr.end(), 1, [](int a, int b) { return a * (2 * b + 1); });
    std::array<double, 3> qstep{};
    std::array<int, 3> nq_arr_in{};
    int ndim = 0;
    for (size_t i = 0; i != 3; ++i)
    {
        if (nq_arr[i] != 0)
        {
            qstep[i] = 1.0 / (2 * nq_arr[i] + 1);
            ndim += 1;
        }
        nq_arr_in[i] = static_cast<int>(nq_arr[i] / a_rate);
    }

    double integ = 0.0;
    int tot_iter = 0;
    for (size_t iter = 0; iter != niter; ++iter)
    {
        double integ_iter = 0.0;
        for (int ig1 = -nq_arr[0]; ig1 != nq_arr[0] + 1; ++ig1)
            for (int ig2 = -nq_arr[1]; ig2 != nq_arr[1] + 1; ++ig2)
                for (int ig3 = -nq_arr[2]; ig3 != nq_arr[2] + 1; ++ig3)
                {
                    if (std::abs(ig1) <= nq_arr_in[0] && std::abs(ig2) <= nq_arr_in[1] && std::abs(ig3) <= nq_arr_in[2])
                        continue;
                    ModuleBase::Vector3<double> qvec; // direct
                    qvec.x = qstep[0] * ig1;
                    qvec.y = qstep[1] * ig2;
                    qvec.z = qstep[2] * ig3;
                    integ_iter += func_cal_fq(qvec * G);
                }
        integ_iter /= nqs * pow(a_rate, ndim * iter); // Each iteration reduces dq by a
                                                      // multiple of a_rate
        integ += integ_iter;
        if (iter != 0 && integ_iter < eps)
            break;
        std::for_each(qstep.begin(), qstep.end(), [&a_rate](double& qs) { qs /= a_rate; });
        tot_iter += 1;
    }

    if (tot_iter == niter)
        ModuleBase::WARNING_QUIT("Singular_Value::Iter_Integral", "not converged!");

    return integ;
}
} // namespace Singular_Value

#endif