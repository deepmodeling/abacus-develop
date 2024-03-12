//==========================================================
// AUTHOR : Yuyang Ji
// DATE : 2019-04-22
// UPDATE : 2021-4-19
//==========================================================

#include "vdwd3_parameters.h"

#include "module_base/constants.h"

namespace vdw
{

void Vdwd3Parameters::initial_parameters(const Input &input)
{
    mxc_.resize(max_elem_, 1);
    r0ab_.resize(max_elem_, std::vector<double>(max_elem_, 0.0));

    c6ab_.resize(3,
                 std::vector<std::vector<std::vector<std::vector<double>>>>(
                     5,
                     std::vector<std::vector<std::vector<double>>>(
                         5,
                         std::vector<std::vector<double>>(max_elem_, std::vector<double>(max_elem_, 0.0)))));

    s6_ = std::stod(input.vdw_s6);
    s18_ = std::stod(input.vdw_s8);
    rs6_ = std::stod(input.vdw_a1);
    rs18_ = std::stod(input.vdw_a2);
    abc_ = input.vdw_abc;
    version_ = input.vdw_method;
    model_ = input.vdw_cutoff_type;
    if (input.vdw_cutoff_type == "radius")
    {
        if (input.vdw_radius_unit == "Bohr")
        {
            rthr2_ = std::pow(std::stod(input.vdw_cutoff_radius), 2);
        }
        else
        {
            rthr2_ = std::pow((std::stod(input.vdw_cutoff_radius) / ModuleBase::BOHR_TO_A), 2);
        }
        if (input.vdw_cn_thr_unit == "Bohr")
        {
            cn_thr2_ = std::pow(input.vdw_cn_thr, 2);
        }
        else
        {
            cn_thr2_ = std::pow((input.vdw_cn_thr / ModuleBase::BOHR_TO_A), 2);
        }
    }
    else if (input.vdw_cutoff_type == "period")
    {
        period_ = input.vdw_cutoff_period;
    }
    init_C6();
    init_r2r4();
    init_rcov();
    init_r0ab();
}

int Vdwd3Parameters::limit(int &i)
{
    int icn = 1;
    while (i >= 100)
    {
        i -= 100;
        icn += 1;
    }
    return icn;
}

} // namespace vdw