/**
 * @file esolver_dp.cpp
 * @brief Implementation of ESolver_DP class for DeePMD method.
 *
 * This file contains the implementation of the ESolver_DP class, which is used for solving the energy and forces in a
 * Deep Potential Molecular Dynamics (DeePMD) simulation.
 * DeePMD is a method for training deep neural networks to accurately predict the potential energy surface of a
 * molecular system.
 *
 * For more information about DeePMD, see the following reference:
 *
 * Han Wang, Linfeng Zhang, Jiequn Han, and Roberto Car.
 * "DeePMD-kit: A deep learning package for many-body potential energy representation and molecular dynamics,"
 * Computer Physics Communications 228, 178-184 (2018). https://doi.org/10.1016/j.cpc.2018.03.016
 *
 * @author YuLiu98
 * @date 2023-05-15
 */
#include "esolver_dp.h"

#include "module_base/parallel_common.h"
#include "module_base/timer.h"
#include "module_io/output_log.h"

#include <iomanip>
#include <sstream>
#include <unordered_map>

namespace ModuleESolver
{

void ESolver_DP::before_all_runners(const Input_para& inp, UnitCell& ucell)
{
    ucell_ = &ucell;
    dp_potential = 0;
    dp_force.create(ucell.nat, 3);
    dp_virial.create(3, 3);

    cell.resize(9);
    atype.resize(ucell.nat);
    coord.resize(3 * ucell.nat);

    fparam = inp.mdp.dp_fparam;
    aparam = inp.mdp.dp_aparam;

#ifdef __DPMD
    /// determine the type map from STRU to DP model
    type_map(ucell);
#endif
}

void ESolver_DP::runner(const int istep, UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_DP", "runner");
    ModuleBase::timer::tick("ESolver_DP", "runner");

    cell[0] = ucell.latvec.e11 * ucell.lat0_angstrom;
    cell[1] = ucell.latvec.e12 * ucell.lat0_angstrom;
    cell[2] = ucell.latvec.e13 * ucell.lat0_angstrom;
    cell[3] = ucell.latvec.e21 * ucell.lat0_angstrom;
    cell[4] = ucell.latvec.e22 * ucell.lat0_angstrom;
    cell[5] = ucell.latvec.e23 * ucell.lat0_angstrom;
    cell[6] = ucell.latvec.e31 * ucell.lat0_angstrom;
    cell[7] = ucell.latvec.e32 * ucell.lat0_angstrom;
    cell[8] = ucell.latvec.e33 * ucell.lat0_angstrom;

    int iat = 0;
    for (int it = 0; it < ucell.ntype; ++it)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
        {
            coord[3 * iat] = ucell.atoms[it].tau[ia].x * ucell.lat0_angstrom;
            coord[3 * iat + 1] = ucell.atoms[it].tau[ia].y * ucell.lat0_angstrom;
            coord[3 * iat + 2] = ucell.atoms[it].tau[ia].z * ucell.lat0_angstrom;
            iat++;
        }
    }
    assert(ucell.nat == iat);

#ifdef __DPMD
    std::vector<double> f, v;
    dp_potential = 0;
    dp_force.zero_out();
    dp_virial.zero_out();

    dp.compute(dp_potential, f, v, coord, atype, cell, fparam, aparam);

    dp_potential /= ModuleBase::Ry_to_eV;
    GlobalV::ofs_running << " final etot is " << std::setprecision(11) << dp_potential * ModuleBase::Ry_to_eV << " eV"
                         << std::endl;

    const double fact_f = ModuleBase::Ry_to_eV * ModuleBase::ANGSTROM_AU;
    const double fact_v = ucell.omega * ModuleBase::Ry_to_eV;

    for (int i = 0; i < ucell.nat; ++i)
    {
        dp_force(i, 0) = f[3 * i] / fact_f;
        dp_force(i, 1) = f[3 * i + 1] / fact_f;
        dp_force(i, 2) = f[3 * i + 2] / fact_f;
    }

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            dp_virial(i, j) = v[3 * i + j] / fact_v;
        }
    }
#else
    ModuleBase::WARNING_QUIT("ESolver_DP", "Please recompile with -D__DPMD");
#endif
    ModuleBase::timer::tick("ESolver_DP", "runner");
}

double ESolver_DP::cal_energy()
{
    return dp_potential;
}

void ESolver_DP::cal_force(ModuleBase::matrix& force)
{
    force = dp_force;
    ModuleIO::print_force(GlobalV::ofs_running, *ucell_, "TOTAL-FORCE (eV/Angstrom)", force, false);
}

void ESolver_DP::cal_stress(ModuleBase::matrix& stress)
{
    stress = dp_virial;

    // external stress
    double unit_transform = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
    double external_stress[3] = {PARAM.inp.press1, PARAM.inp.press2, PARAM.inp.press3};
    for (int i = 0; i < 3; i++)
    {
        stress(i, i) -= external_stress[i] / unit_transform;
    }

    ModuleIO::print_stress("TOTAL-STRESS", stress, true, false);
}

void ESolver_DP::after_all_runners()
{
    GlobalV::ofs_running << "\n\n --------------------------------------------" << std::endl;
    GlobalV::ofs_running << std::setprecision(16);
    GlobalV::ofs_running << " !FINAL_ETOT_IS " << dp_potential * ModuleBase::Ry_to_eV << " eV" << std::endl;
    GlobalV::ofs_running << " --------------------------------------------\n\n" << std::endl;
}

#ifdef __DPMD
void ESolver_DP::type_map(const UnitCell& ucell)
{
    std::string type = "";
    dp.get_type_map(type);
    std::stringstream ss(type);
    std::unordered_map<std::string, int> label;
    std::string temp;
    int index = 0;
    while (ss >> temp)
    {
        label[temp] = index;
        index++;
    }

    std::cout << "\n type map of model file " << dp_file << " " << std::endl;
    std::cout << " ----------------------------------------------------------------";
    int count = 0;
    for (auto it = label.begin(); it != label.end(); ++it)
    {
        if (count % 5 == 0)
        {
            std::cout << std::endl;
            std::cout << "  ";
        }
        count++;
        temp = it->first + ": " + std::to_string(it->second);
        std::cout << std::left << std::setw(10) << temp;
    }
    std::cout << "\n -----------------------------------------------------------------" << std::endl;

    int iat = 0;
    for (int it = 0; it < ucell.ntype; ++it)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
        {
            if (label.find(ucell.atoms[it].label) == label.end())
            {
                ModuleBase::WARNING_QUIT("ESolver_DP",
                                         "The label " + ucell.atoms[it].label + " is not found in the type map.");
            }
            atype[iat] = label[ucell.atoms[it].label];
            // if (ia == 0)
            //     std::cout << "type: " << atype[iat] << std::endl;
            iat++;
        }
    }
    assert(ucell.nat == iat);
}
#endif
} // namespace ModuleESolver
