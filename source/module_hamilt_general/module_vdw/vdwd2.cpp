//==========================================================
// AUTHOR : Peize Lin
// DATE : 2014-04-25
// UPDATE : 2019-04-26
//==========================================================

#include "vdwd2.h"
#include "module_base/timer.h"

namespace vdw
{

void Vdwd2::cal_energy()
{
    ModuleBase::TITLE("Vdwd2", "energy");
    ModuleBase::timer::tick("Vdwd2", "energy");
    para_.initset(ucell_);
    energy_ = 0;

    auto energy = [&](double r,
                      double R0_sum,
                      double C6_product,
                      double r_sqr,
                      int,
                      int,
                      const ModuleBase::Vector3<double> &,
                      const ModuleBase::Vector3<double> &) {
        const double tmp_damp_recip = 1 + exp(-para_.damping() * (r / R0_sum - 1));
        energy_ -= C6_product / pow(r_sqr, 3) / tmp_damp_recip / 2;
    };
    index_loops(energy);
    energy_ *= para_.scaling();
    ModuleBase::timer::tick("Vdwd2", "energy");
}

void Vdwd2::cal_force()
{
    ModuleBase::TITLE("Vdwd2", "force");
    ModuleBase::timer::tick("Vdwd2", "force");
    para_.initset(ucell_);
    force_.clear();
    force_.resize(ucell_.nat);

    auto force = [&](double r,
                     double R0_sum,
                     double C6_product,
                     double r_sqr,
                     int it1,
                     int ia1,
                     const ModuleBase::Vector3<double> &tau1,
                     const ModuleBase::Vector3<double> &tau2) {
        const double tmp_exp = exp(-para_.damping() * (r / R0_sum - 1));
        const double tmp_factor = C6_product / pow(r_sqr, 3) / r / (1 + tmp_exp)
                                  * (-6 / r + tmp_exp / (1 + tmp_exp) * para_.damping() / R0_sum);
        force_[ucell_.itia2iat(it1, ia1)] += tmp_factor * (tau1 - tau2);
    };

    index_loops(force);
    std::for_each(force_.begin(), force_.end(), [&](ModuleBase::Vector3<double> &f) {
        f *= para_.scaling() / ucell_.lat0;
    });
    ModuleBase::timer::tick("Vdwd2", "force");
}

void Vdwd2::cal_stress()
{
    ModuleBase::TITLE("Vdwd2", "stress");
    ModuleBase::timer::tick("Vdwd2", "stress");
    para_.initset(ucell_);
    stress_.Zero();

    auto stress = [&](double r,
                      double R0_sum,
                      double C6_product,
                      double r_sqr,
                      int it1,
                      int ia1,
                      const ModuleBase::Vector3<double> &tau1,
                      const ModuleBase::Vector3<double> &tau2) {
        const double tmp_exp = exp(-para_.damping() * (r / R0_sum - 1));
        const double tmp_factor = C6_product / pow(r_sqr, 3) / r / (1 + tmp_exp)
                                  * (-6 / r + tmp_exp / (1 + tmp_exp) * para_.damping() / R0_sum);
        const ModuleBase::Vector3<double> dr = tau2 - tau1;
        stress_ += tmp_factor / 2
                   * ModuleBase::Matrix3(dr.x * dr.x, dr.x * dr.y, dr.x * dr.z,
                                         dr.y * dr.x, dr.y * dr.y, dr.y * dr.z,
                                         dr.z * dr.x, dr.z * dr.y, dr.z * dr.z);
    };

    index_loops(stress);
    stress_ *= para_.scaling() / ucell_.omega;
    ModuleBase::timer::tick("Vdwd2", "stress");
}

} // namespace vdw
