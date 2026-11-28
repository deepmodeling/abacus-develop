#include "dfpt_symm.h"

#include "source_base/constants.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"

namespace ModuleDFPT
{

std::vector<int> little_group_of_q(
    const ModuleBase::Matrix3 &latvec,
    const ModuleBase::Matrix3 *gmatrix,
    int nops,
    const ModuleBase::Vector3<double> &q_cart)
{
    ModuleBase::TITLE("DFPT_Symm", "little_group_of_q");
    ModuleBase::timer::start("DFPT_Symm", "little_group_of_q");
    std::vector<int> group;
    if (gmatrix == nullptr || nops <= 0)
    {
        ModuleBase::timer::end("DFPT_Symm", "little_group_of_q");
        return group;
    }
    // cartesian rotation images of the fractional operations, following
    // the symmetrize_mat3 construction (R = latvec^-1 G latvec)
    const ModuleBase::Matrix3 inv_lat = latvec.Inverse();
    // direct primitive vectors (rows of latvec) for the reciprocal test
    const ModuleBase::Vector3<double> a1(latvec.e11, latvec.e12, latvec.e13);
    const ModuleBase::Vector3<double> a2(latvec.e21, latvec.e22, latvec.e23);
    const ModuleBase::Vector3<double> a3(latvec.e31, latvec.e32, latvec.e33);
    const double tol = 1.0e-6; ///< empirical parameter: integrality roundoff of the double dot products
    for (int iop = 0; iop < nops; ++iop)
    {
        const ModuleBase::Matrix3 rot = inv_lat * gmatrix[iop] * latvec;
        const ModuleBase::Vector3<double> rq = rot * q_cart;
        for (int s = 1; s >= -1; s -= 2)
        {
            const ModuleBase::Vector3<double> diff = (s == 1) ? (rq - q_cart) : (rq + q_cart);
            const double d1 = diff * a1 / ModuleBase::TWO_PI;
            const double d2 = diff * a2 / ModuleBase::TWO_PI;
            const double d3 = diff * a3 / ModuleBase::TWO_PI;
            const bool recip = std::abs(d1 - std::round(d1)) < tol && std::abs(d2 - std::round(d2)) < tol &&
                               std::abs(d3 - std::round(d3)) < tol;
            if (recip)
            {
                group.push_back(iop);
                break; // either branch suffices; keep each operation once
            }
        }
    }
    ModuleBase::timer::end("DFPT_Symm", "little_group_of_q");
    return group;
}

} // namespace ModuleDFPT
