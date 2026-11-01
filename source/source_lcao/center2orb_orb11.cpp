//=========================================================
// AUTHOR : Peize Lin
// DATE : 2016-01-24
//=========================================================

#include "center2orb_orb11.h"

#include "center2orb.h"
#include "source_base/constants.h"
#include "source_base/math_polyint.h"
#include "source_base/sph_bessel_recursive.h"
#include "source_base/ylm.h"

#include <cmath>

Center2_Orb::Orb11::Orb11(const Numerical_Orbital_Lm& nA_in,
                          const Numerical_Orbital_Lm& nB_in,
                          const ModuleBase::Sph_Bessel_Recursive::D2* psb,
                          const ORB_gaunt_table& MGT_in)
    : nA(nA_in), nB(nB_in), psb_(psb), MGT(MGT_in)
{
}

void Center2_Orb::Orb11::init_radial_table()
{
    const int rmesh = Center2_Orb::get_rmesh(this->nA.getRcut(), this->nB.getRcut(), dr_);

    const int LA = this->nA.getL();
    const int LB = this->nB.getL();
    for (int LAB = std::abs(LA - LB); LAB <= LA + LB; ++LAB)
    {
        if ((LAB - std::abs(LA - LB)) % 2 == 1) // if LA+LB-LAB == odd, then Gaunt_Coefficients = 0
        {
            continue;
        }

        this->Table_r[LAB].resize(rmesh, 0);
        this->Table_dr[LAB].resize(rmesh, 0);

        Center2_Orb::cal_ST_Phi12_R(1,
                                    LAB,
                                    this->nA,
                                    this->nB,
                                    rmesh,
                                    this->Table_r[LAB],
                                    this->Table_dr[LAB],
                                    psb_);
    }
    return;
}

void Center2_Orb::Orb11::init_radial_table(const std::set<size_t>& radials)
{
    const int LA = this->nA.getL();
    const int LB = this->nB.getL();

    const size_t rmesh = Center2_Orb::get_rmesh(this->nA.getRcut(), this->nB.getRcut(), dr_);

    std::set<size_t> radials_used;
    for (const size_t& ir: radials)
    {
        if (ir < rmesh)
        {
            radials_used.insert(ir);
        }
    }

    for (int LAB = std::abs(LA - LB); LAB <= LA + LB; ++LAB)
    {
        if ((LAB - std::abs(LA - LB)) % 2 == 1) // if LA+LB-LAB == odd, then Gaunt_Coefficients = 0
        { 
            continue;
        }

        this->Table_r[LAB].resize(rmesh, 0);
        this->Table_dr[LAB].resize(rmesh, 0);

        Center2_Orb::cal_ST_Phi12_R(1,
                                    LAB,
                                    this->nA,
                                    this->nB,
                                    radials_used,
                                    this->Table_r[LAB],
                                    this->Table_dr[LAB],
                                    psb_);
    }
}

double Center2_Orb::Orb11::cal_overlap(const ModuleBase::Vector3<double>& RA,
                                       const ModuleBase::Vector3<double>& RB,
                                       const int& mA,
                                       const int& mB) const
{
    const double tiny1 = 1e-12; // why 1e-12?
    const double tiny2 = 1e-10; // why 1e-10?

    const ModuleBase::Vector3<double> delta_R = RB - RA;
    const double distance_true = delta_R.norm();
    const double distance = (distance_true >= tiny1) ? distance_true : distance_true + tiny1;
    const double RcutA = this->nA.getRcut();
    const double RcutB = this->nB.getRcut();
    if (distance > (RcutA + RcutB))
    {
        return 0.0;
    }

    const int LA = this->nA.getL();
    const int LB = this->nB.getL();

    std::vector<double> rly;
    ModuleBase::Ylm::rl_sph_harm(LA + LB, // max LAB
                                 delta_R.x,
                                 delta_R.y,
                                 delta_R.z,
                                 rly);

    double overlap = 0.0;
    const int idx1 = this->MGT.get_lm_index(LA, mA);
    const int idx2 = this->MGT.get_lm_index(LB, mB);
    const double* Gaunt_Coefficients_ptr = &(this->MGT.Gaunt_Coefficients(idx1, idx2, 0));

    const int rmesh = Center2_Orb::get_rmesh(RcutA, RcutB, dr_);
    for (const auto& tb_r: this->Table_r)
    {
        const int LAB = tb_r.first;
        const double distance_power = std::pow(distance, LAB);
        const int phase = (LA - LB - LAB) / 2;
        const double i_exp = std::pow(-1.0, phase);
        const double Interp_Tlm = (distance > tiny2) ? ModuleBase::PolyInt::Polynomial_Interpolation(
                                                           tb_r.second.data(),
                                                           rmesh,
                                                           dr_,
                                                           distance)
                                                     : tb_r.second.at(0);

        for (int mAB = 0; mAB != 2 * LAB + 1; ++mAB)
        // const int mAB = mA + mB;
        {
            const int idx3 = this->MGT.get_lm_index(LAB, mAB);
            const double Gaunt_real_A_B_AB = *(Gaunt_Coefficients_ptr + idx3);
            if (0 == Gaunt_real_A_B_AB)
            {
                continue;
            }

            const double ylm_solid = rly[idx3];
            if (0 == ylm_solid)
            {
                continue;
            }
            const double ylm_real = (distance > tiny2) ? ylm_solid / distance_power : ylm_solid;

            overlap +=
                //			pow(2*PI,1.5)
                i_exp * Gaunt_real_A_B_AB * Interp_Tlm * ylm_real;
        }
    }

    return overlap;
}

ModuleBase::Vector3<double> Center2_Orb::Orb11::cal_grad_overlap( // caoyu add 2021-11-19
    const ModuleBase::Vector3<double>& RA,
    const ModuleBase::Vector3<double>& RB,
    const int& mA,
    const int& mB) const
{
    const double tiny1 = 1e-12; // same as `cal_overlap`
    const double tiny2 = 1e-10; // same as `cal_overlap`

    const ModuleBase::Vector3<double> delta_R = RB - RA;
    const double distance_true = delta_R.norm();
    const double distance = (distance_true >= tiny1) ? distance_true : distance_true + tiny1;
    const double RcutA = this->nA.getRcut();
    const double RcutB = this->nB.getRcut();
    if (distance > (RcutA + RcutB))
    {
        return ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
    }

    const int LA = this->nA.getL();
    const int LB = this->nB.getL();
    const int idx1 = this->MGT.get_lm_index(LA, mA);
    const int idx2 = this->MGT.get_lm_index(LB, mB);
    const double* Gaunt_Coefficients_ptr = &(this->MGT.Gaunt_Coefficients(idx1, idx2, 0));

    const int LAB2 = (LA + LB + 1) * (LA + LB + 1);
    std::vector<double> rly(LAB2);
    std::vector<ModuleBase::Vector3<double>> grly;
    std::vector<double> tmp_grly(LAB2 * 3);
    ModuleBase::Ylm::grad_rl_sph_harm(LA + LB, delta_R.x, delta_R.y, delta_R.z, rly.data(), tmp_grly.data());
    grly.reserve(LAB2);
    for (int i=0; i<LAB2; ++i)
    {
        grly.emplace_back(tmp_grly[i*3], tmp_grly[i*3+1], tmp_grly[i*3+2]);
    }

    ModuleBase::Vector3<double> grad_overlap(0.0, 0.0, 0.0);

    const int rmesh = Center2_Orb::get_rmesh(RcutA, RcutB, dr_);
    for (const auto& tb_r: this->Table_r)
    {
        const int LAB = tb_r.first;
        const double distance_power = std::pow(distance, LAB);
        const int phase = (LA - LB - LAB) / 2;
        const double i_exp = std::pow(-1.0, phase);
        const double Interp_Tlm = (distance > tiny2) ? ModuleBase::PolyInt::Polynomial_Interpolation(
                                                           tb_r.second.data(),
                                                           rmesh,
                                                           dr_,
                                                           distance)
                                                     : tb_r.second.at(0);
        const double grad_Interp_Tlm = (distance > tiny2) ? ModuleBase::PolyInt::Polynomial_Interpolation(
                                                                this->Table_dr.at(LAB).data(),
                                                                rmesh,
                                                                dr_,
                                                                distance) -
                                                                Interp_Tlm * LAB / distance
                                                          : 0.0;
        for (int mAB = 0; mAB != 2 * LAB + 1; ++mAB)
        // const int mAB = mA + mB;
        {
            const int idx3 = this->MGT.get_lm_index(LAB, mAB);
            const double Gaunt_real_A_B_AB = *(Gaunt_Coefficients_ptr + idx3);
            if (0 == Gaunt_real_A_B_AB)
            {
                continue;
            }

            const double ylm_solid = rly[idx3];
            const double ylm_real = (distance > tiny2) ? ylm_solid / distance_power : ylm_solid;

            const ModuleBase::Vector3<double> gylm_solid = grly[idx3];
            const ModuleBase::Vector3<double> gylm_real = (distance > tiny2) ? gylm_solid / distance_power : gylm_solid;

            grad_overlap += i_exp //			pow(2*PI,1.5)
                            * Gaunt_real_A_B_AB
                            * (Interp_Tlm * gylm_real + grad_Interp_Tlm * ylm_real * delta_R / distance);
        }
    }

    return grad_overlap;
}

template <typename T> std::vector<T> Center2_Orb::Orb11::contract_angular(const std::vector<T> &components) const
{
    const int LA = this->nA.getL();
    const int LB = this->nB.getL();
    const int count_A = 2 * LA + 1;
    const int count_B = 2 * LB + 1;
    const int count = count_A * count_B;
    std::vector<T> result(count);
    for (int ma = 0; ma < count_A; ++ma)
    {
        const int idx_A = this->MGT.get_lm_index(LA, ma);
        for (int mb = 0; mb < count_B; ++mb)
        {
            const int idx_B = this->MGT.get_lm_index(LB, mb);
            const double *coefficients = &this->MGT.Gaunt_Coefficients(idx_A, idx_B, 0);
            T &value = result[ma * count_B + mb];
            // Keep the scalar path's angular-momentum accumulation order.
            for (const std::pair<const int, std::vector<double>> &table : this->Table_r)
            {
                const int l = table.first;
                for (int m = 0; m <= 2 * l; ++m)
                {
                    const int idx = this->MGT.get_lm_index(l, m);
                    if (coefficients[idx] != 0.0)
                    {
                        value += coefficients[idx] * components[idx];
                    }
                }
            }
        }
    }
    return result;
}

std::vector<double> Center2_Orb::Orb11::cal_overlap(
    const ModuleBase::Vector3<double> &RA, const ModuleBase::Vector3<double> &RB) const
{
    const double tiny1 = 1e-12;
    const double tiny2 = 1e-10;
    const ModuleBase::Vector3<double> delta_R = RB - RA;
    const double distance_true = delta_R.norm();
    const double distance = distance_true >= tiny1 ? distance_true : distance_true + tiny1;
    const int LA = this->nA.getL();
    const int LB = this->nB.getL();
    const double RcutA = this->nA.getRcut();
    const double RcutB = this->nB.getRcut();
    if (distance > RcutA + RcutB)
    {
        const int count = (2 * LA + 1) * (2 * LB + 1);
        return std::vector<double>(count);
    }

    const int lmax = LA + LB;
    std::vector<double> rly;
    ModuleBase::Ylm::rl_sph_harm(
        lmax,
        delta_R.x,
        delta_R.y,
        delta_R.z,
        rly);
    const int angular_count = (lmax + 1) * (lmax + 1);
    std::vector<double> components(angular_count);
    const int rmesh = Center2_Orb::get_rmesh(RcutA, RcutB, this->dr_);
    for (const std::pair<const int, std::vector<double>> &table : this->Table_r)
    {
        const int l = table.first;
        const double distance_power = std::pow(distance, l);
        const int phase = (LA - LB - l) / 2;
        const double i_exp = std::pow(-1.0, phase);
        const double radial = distance > tiny2 ? ModuleBase::PolyInt::Polynomial_Interpolation(
                                                     table.second.data(),
                                                     rmesh,
                                                     dr_,
                                                     distance)
                                               : table.second.at(0);
        for (int m = 0; m <= 2 * l; ++m)
        {
            const int idx = this->MGT.get_lm_index(l, m);
            const double angular = distance > tiny2 ? rly[idx] / distance_power : rly[idx];
            components[idx] = i_exp * radial * angular;
        }
    }
    return this->contract_angular(components);
}

std::vector<ModuleBase::Vector3<double>> Center2_Orb::Orb11::cal_grad_overlap(
    const ModuleBase::Vector3<double> &RA, const ModuleBase::Vector3<double> &RB) const
{
    const double tiny1 = 1e-12;
    const double tiny2 = 1e-10;
    const ModuleBase::Vector3<double> delta_R = RB - RA;
    const double distance_true = delta_R.norm();
    const double distance = distance_true >= tiny1 ? distance_true : distance_true + tiny1;
    const int LA = this->nA.getL();
    const int LB = this->nB.getL();
    const double RcutA = this->nA.getRcut();
    const double RcutB = this->nB.getRcut();
    if (distance > RcutA + RcutB)
    {
        const int count = (2 * LA + 1) * (2 * LB + 1);
        return std::vector<ModuleBase::Vector3<double>>(count);
    }

    const int lmax = LA + LB;
    const int angular_count = (lmax + 1) * (lmax + 1);
    const int gradient_count = 3 * angular_count;
    std::vector<double> rly(angular_count);
    std::vector<double> grly(gradient_count);
    ModuleBase::Ylm::grad_rl_sph_harm(
        lmax,
        delta_R.x,
        delta_R.y,
        delta_R.z,
        rly.data(),
        grly.data());
    std::vector<ModuleBase::Vector3<double>> components(angular_count);
    const int rmesh = Center2_Orb::get_rmesh(RcutA, RcutB, this->dr_);
    for (const std::pair<const int, std::vector<double>> &table : this->Table_r)
    {
        const int l = table.first;
        const double distance_power = std::pow(distance, l);
        const int phase = (LA - LB - l) / 2;
        const double i_exp = std::pow(-1.0, phase);
        const double radial = distance > tiny2 ? ModuleBase::PolyInt::Polynomial_Interpolation(
                                                     table.second.data(),
                                                     rmesh,
                                                     dr_,
                                                     distance)
                                               : table.second.at(0);
        const double derivative = distance > tiny2 ? ModuleBase::PolyInt::Polynomial_Interpolation(
                                                         this->Table_dr.at(l).data(),
                                                         rmesh,
                                                         dr_,
                                                         distance) -
                                                         radial * l / distance
                                                   : 0.0;
        for (int m = 0; m <= 2 * l; ++m)
        {
            const int idx = this->MGT.get_lm_index(l, m);
            const double angular = distance > tiny2 ? rly[idx] / distance_power : rly[idx];
            const int offset = 3 * idx;
            const double gx = grly[offset];
            const double gy = grly[offset + 1];
            const double gz = grly[offset + 2];
            const ModuleBase::Vector3<double> solid_gradient(gx, gy, gz);
            const ModuleBase::Vector3<double> angular_gradient =
                distance > tiny2 ? solid_gradient / distance_power : solid_gradient;
            components[idx] = i_exp * (radial * angular_gradient + derivative * angular * delta_R / distance);
        }
    }
    return this->contract_angular(components);
}
