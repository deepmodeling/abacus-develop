//==========================================================
// AUTHOR : Yuyang Ji
// DATE : 2019-04-22
// UPDATE : 2021-4-19
//==========================================================

#include "vdwd3_parameters.h"
#include "module_base/constants.h"
#include <map>
#include "dftd3.h"

namespace vdw
{

void Vdwd3Parameters::initial_parameters(const Input_para &input)
{
    // initialize the dftd3 parameters
    mxc_.resize(max_elem_, 1);
    r0ab_.resize(max_elem_, std::vector<double>(max_elem_, 0.0));

    c6ab_.resize(3,
                 std::vector<std::vector<std::vector<std::vector<double>>>>(
                     5,
                     std::vector<std::vector<std::vector<double>>>(
                         5,
                         std::vector<std::vector<double>>(max_elem_, std::vector<double>(max_elem_, 0.0)))));
    const std::string xc = input.dft_functional;
    const std::string vdw_method = input.vdw_method;
    const std::map<std::string, std::string> dftd3_method = {{"d3_bj", "bj"}, {"d3_0", "zero"},
                                                             {"d3_bjm", "bjm"}, {"d3_0m", "zerom"},
                                                             {"op", "op"}};
    std::vector<double> param;
    std::cout << " VDW: search for DFTD3 parameters for " << xc << " with " << vdw_method << std::endl;
    DFTD3::search(xc, dftd3_method.at(vdw_method), param);
    s6_ = param[0];
    s18_ = param[3];
    rs6_ = param[1];
    rs18_ = param[5];
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