//==========================================================
// AUTHOR : Peize Lin
// DATE : 2014-04-25
// UPDATE : 2019-04-26
//==========================================================

#include "vdwd2.h"
#include "source_base/timer.h"

#include <algorithm>
#include <cmath>

namespace vdw
{

void Vdwd2::evaluate_impl(const VdwRequest& request, VdwResult& result)
{
    ModuleBase::TITLE("Vdwd2", "evaluate");
    ModuleBase::timer::start("Vdwd2", "evaluate");

    para_.initset(ucell_);

    if (request.force)
    {
        result.force.resize(ucell_.nat);
    }
    if (request.stress)
    {
        result.stress.Zero();
    }

    const bool need_derivatives = request.force || request.stress;
    auto evaluate_pair = [&](double r,
                             double R0_sum,
                             double C6_product,
                             double r_sqr,
                             int it1,
                             int ia1,
                             const ModuleBase::Vector3<double>& tau1,
                             const ModuleBase::Vector3<double>& tau2) {
        const double tmp_exp = exp(-para_.damping() * (r / R0_sum - 1));
        const double tmp_damp_recip = 1.0 + tmp_exp;

        result.energy -= C6_product / pow(r_sqr, 3) / tmp_damp_recip / 2.0;

        if (!need_derivatives)
        {
            return;
        }

        const double tmp_factor = C6_product / pow(r_sqr, 3) / r / tmp_damp_recip
                                  * (-6.0 / r
                                     + tmp_exp / tmp_damp_recip * para_.damping() / R0_sum);

        if (request.force)
        {
            result.force[ucell_.itia2iat(it1, ia1)] += tmp_factor * (tau1 - tau2);
        }

        if (request.stress)
        {
            const ModuleBase::Vector3<double> dr = tau2 - tau1;
            result.stress += tmp_factor / 2.0
                             * ModuleBase::Matrix3(dr.x * dr.x,
                                                   dr.x * dr.y,
                                                   dr.x * dr.z,
                                                   dr.y * dr.x,
                                                   dr.y * dr.y,
                                                   dr.y * dr.z,
                                                   dr.z * dr.x,
                                                   dr.z * dr.y,
                                                   dr.z * dr.z);
        }
    };

    index_loops(evaluate_pair);

    result.energy *= para_.scaling();

    if (request.force)
    {
        std::for_each(result.force.begin(), result.force.end(), [&](ModuleBase::Vector3<double>& force) {
            force *= para_.scaling() / ucell_.lat0;
        });
        result.has_force = true;
    }

    if (request.stress)
    {
        result.stress *= para_.scaling() / ucell_.omega;
        result.has_stress = true;
    }

    ModuleBase::timer::end("Vdwd2", "evaluate");
}

} // namespace vdw
