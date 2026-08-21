#include "vdwd3.h"

#include "vdwd3_evaluator.h"
#include "vdwd3_parameters.h"
#include "source_base/constants.h"
#include "source_base/element_name.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <iomanip>
#include <string>

namespace vdw
{
namespace
{

int atomic_number_from_symbol(const std::string& symbol)
{
    for (std::size_t index = 0; index < ModuleBase::element_name.size(); ++index)
    {
        if (symbol == ModuleBase::element_name[index])
        {
            return static_cast<int>(index) + 1;
        }
    }
    ModuleBase::WARNING_QUIT("Vdwd3::atomic_number_from_symbol", "Unknown element symbol: " + symbol);
}

double length_to_bohr(double value, const std::string& unit)
{
    if (unit == "Bohr")
    {
        return value;
    }
    if (unit == "A")
    {
        return value / ModuleBase::BOHR_TO_A;
    }
    ModuleBase::WARNING_QUIT("Vdwd3::length_to_bohr", "Unsupported length unit: " + unit);
}

double cutoff_to_bohr(const std::string& value, const std::string& unit)
{
    return length_to_bohr(std::stod(value), unit);
}

d3::Vec3 to_d3_vector(const ModuleBase::Vector3<double>& value)
{
    return d3::Vec3(value.x, value.y, value.z);
}

const char* damping_name(d3::Damping damping)
{
    return damping == d3::Damping::Rational ? "rational (BJ)" : "zero";
}

} // namespace

Vdwd3::Vdwd3(const UnitCell& unit_in,
             const std::string& xc_name,
             const Input_para& input,
             std::ofstream* plog)
    : Vdw(unit_in)
{
    if (input.vdw_cutoff_type != "radius")
    {
        ModuleBase::WARNING_QUIT("Vdwd3::Vdwd3",
                                 "DFT-D3 requires vdw_cutoff_type=radius; "
                                 "the legacy period-count cutoff is no longer supported.");
    }

    if (input.vdw_method == "d3_bj")
    {
        parameters_.damping = d3::Damping::Rational;
    }
    else if (input.vdw_method == "d3_0")
    {
        parameters_.damping = d3::Damping::Zero;
    }
    else
    {
        ModuleBase::WARNING_QUIT("Vdwd3::Vdwd3", "Unsupported DFT-D3 method: " + input.vdw_method);
    }

    const bool all_custom = input.vdw_s6 != "default"
                            && input.vdw_s8 != "default"
                            && input.vdw_a1 != "default"
                            && input.vdw_a2 != "default";
    if (all_custom)
    {
        canonical_method_ = d3::canonicalize_method_name(xc_name);
    }
    else if (!d3::lookup_parameters(xc_name,
                                    parameters_.damping,
                                    parameters_,
                                    canonical_method_))
    {
        ModuleBase::WARNING_QUIT("Vdwd3::Vdwd3",
                                 "No s-dftd3 damping parameters found for XC functional '"
                                     + xc_name + "'. Define vdw_s6, vdw_s8, vdw_a1 and vdw_a2 "
                                                 "together to use a fully custom parameter set.");
    }

    if (input.vdw_s6 != "default")
    {
        parameters_.s6 = std::stod(input.vdw_s6);
    }
    if (input.vdw_s8 != "default")
    {
        parameters_.s8 = std::stod(input.vdw_s8);
    }
    if (input.vdw_a1 != "default")
    {
        if (parameters_.damping == d3::Damping::Rational)
        {
            parameters_.a1 = std::stod(input.vdw_a1);
        }
        else
        {
            parameters_.rs6 = std::stod(input.vdw_a1);
        }
    }
    if (input.vdw_a2 != "default")
    {
        if (parameters_.damping == d3::Damping::Rational)
        {
            parameters_.a2 = std::stod(input.vdw_a2);
        }
        else
        {
            parameters_.rs8 = std::stod(input.vdw_a2);
        }
    }
    parameters_.s9 = input.vdw_abc ? 1.0 : 0.0;

    cutoffs_.disp2 = cutoff_to_bohr(input.vdw_cutoff_radius, input.vdw_radius_unit);
    cutoffs_.disp3 = std::min(40.0, cutoffs_.disp2);
    cutoffs_.cn = length_to_bohr(input.vdw_cn_thr, input.vdw_cn_thr_unit);
    cutoffs_.width2 = input.vdw_cutoff_width2;
    cutoffs_.width3 = input.vdw_cutoff_width3;

    if (cutoffs_.width2 < 0.0 || cutoffs_.width2 > cutoffs_.disp2)
    {
        ModuleBase::WARNING_QUIT("Vdwd3::Vdwd3",
                                 "vdw_cutoff_width2 must satisfy "
                                 "0 <= width <= two-body cutoff");
    }
    if (cutoffs_.width3 < 0.0 || cutoffs_.width3 > cutoffs_.disp3)
    {
        ModuleBase::WARNING_QUIT("Vdwd3::Vdwd3",
                                 "vdw_cutoff_width3 must satisfy "
                                 "0 <= width <= three-body cutoff");
    }

    write_parameters(plog);
}

void Vdwd3::write_parameters(std::ofstream* plog) const
{
    if (plog == nullptr)
    {
        return;
    }

    *plog << "\nDFT-D3 parameters (s-dftd3 v1.5.0 numerical specification)\n"
          << "XC functional: " << canonical_method_ << '\n'
          << "damping: " << damping_name(parameters_.damping) << '\n'
          << std::setprecision(10)
          << "s6=" << parameters_.s6 << " s8=" << parameters_.s8
          << " s9=" << parameters_.s9;
    if (parameters_.damping == d3::Damping::Rational)
    {
        *plog << " a1=" << parameters_.a1 << " a2=" << parameters_.a2;
    }
    else
    {
        *plog << " rs6=" << parameters_.rs6 << " rs8=" << parameters_.rs8;
    }
    *plog << "\ncutoffs (Bohr): disp2=" << cutoffs_.disp2
          << " disp3=" << cutoffs_.disp3 << " CN=" << cutoffs_.cn
          << " smooth_width_2b=" << cutoffs_.width2
          << " smooth_width_3b=" << cutoffs_.width3 << std::endl;
}

d3::Structure Vdwd3::build_structure() const
{
    d3::Structure structure;
    structure.atomic_numbers.reserve(ucell_.nat);
    structure.positions.reserve(ucell_.nat);

    for (int it = 0; it < ucell_.ntype; ++it)
    {
        const int atomic_number = atomic_number_from_symbol(ucell_.atoms[it].ncpp.psd);
        for (int ia = 0; ia < ucell_.atoms[it].na; ++ia)
        {
            structure.atomic_numbers.push_back(atomic_number);
            structure.positions.push_back(to_d3_vector(ucell_.atoms[it].tau[ia] * ucell_.lat0));
        }
    }

    structure.lattice = {{to_d3_vector(ucell_.a1 * ucell_.lat0),
                          to_d3_vector(ucell_.a2 * ucell_.lat0),
                          to_d3_vector(ucell_.a3 * ucell_.lat0)}};
    structure.periodic = {{true, true, true}};
    return structure;
}

void Vdwd3::evaluate_impl(const VdwRequest& request, VdwResult& result)
{
    ModuleBase::TITLE("Vdwd3", "evaluate");
    ModuleBase::timer::start("Vdwd3", "evaluate");

    d3::Result d3_result;
    std::string error;
    const bool derivatives = request.force || request.stress;
    if (!d3::evaluate(build_structure(), parameters_, cutoffs_, derivatives, d3_result, error))
    {
        ModuleBase::WARNING_QUIT("Vdwd3::evaluate", error);
    }

    // The native core follows s-dftd3 and returns Hartree-based quantities.
    result.energy = 2.0 * d3_result.energy;

    if (request.force)
    {
        result.force.resize(ucell_.nat);
        for (int iat = 0; iat < ucell_.nat; ++iat)
        {
            result.force[iat] = ModuleBase::Vector3<double>(-2.0 * d3_result.gradient[iat].x,
                                                            -2.0 * d3_result.gradient[iat].y,
                                                            -2.0 * d3_result.gradient[iat].z);
        }
        result.has_force = true;
    }

    if (request.stress)
    {
        const d3::Matrix3& sigma = d3_result.virial;
        result.stress = ModuleBase::Matrix3(2.0 * sigma.value[0][0],
                                            2.0 * sigma.value[0][1],
                                            2.0 * sigma.value[0][2],
                                            2.0 * sigma.value[1][0],
                                            2.0 * sigma.value[1][1],
                                            2.0 * sigma.value[1][2],
                                            2.0 * sigma.value[2][0],
                                            2.0 * sigma.value[2][1],
                                            2.0 * sigma.value[2][2])
                        / ucell_.omega;
        result.has_stress = true;
    }

    ModuleBase::timer::end("Vdwd3", "evaluate");
}

} // namespace vdw
