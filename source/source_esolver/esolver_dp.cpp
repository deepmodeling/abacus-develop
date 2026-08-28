/**
 * @file esolver_dp.cpp
#include "source_io/module_parameter/parameter.h"
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
#include "source_base/parallel_common.h"
#include "source_base/timer.h"
#include "source_io/module_output/output_log.h"
#include "source_io/module_parameter/parameter.h"

#include <iomanip>
#include <sstream>
#include <unordered_map>

using namespace ModuleESolver;

void ESolver_DP::before_all_runners(BaseCell& basecell, const Input_para& inp)
{
    basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    this->inp_ = &inp;

    dp_potential = 0;
    dp_force.create(ucell.nat, 3);
    dp_virial.create(3, 3);
    dp_cell.resize(9);
    dp_coord.resize(3 * ucell.nat);
    dp_model_force.clear();
    dp_model_virial.clear();

    atype.resize(ucell.nat);

    // Build flat atom index for OpenMP coordinate fill in runner()
    atom_type_index.resize(ucell.nat);
    atom_local_index.resize(ucell.nat);
    int iat = 0;
    for (int it = 0; it < ucell.ntype; ++it)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
        {
            atom_type_index[iat] = it;
            atom_local_index[iat] = ia;
            iat++;
        }
    }

    rescaling = inp.mdp.dp_rescaling;
    fparam = inp.mdp.dp_fparam;
    aparam = inp.mdp.dp_aparam;

#ifdef __DPMD
    /// determine the type map from STRU to DP model
    type_map(ucell);
#endif
}

void ESolver_DP::runner(BaseCell& basecell, const int istep)
{
    basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    ModuleBase::TITLE("ESolver_DP", "runner");
    ModuleBase::timer::start("ESolver_DP", "runner");

    dp_cell[0] = ucell.latvec.e11 * ucell.lat0_angstrom;
    dp_cell[1] = ucell.latvec.e12 * ucell.lat0_angstrom;
    dp_cell[2] = ucell.latvec.e13 * ucell.lat0_angstrom;
    dp_cell[3] = ucell.latvec.e21 * ucell.lat0_angstrom;
    dp_cell[4] = ucell.latvec.e22 * ucell.lat0_angstrom;
    dp_cell[5] = ucell.latvec.e23 * ucell.lat0_angstrom;
    dp_cell[6] = ucell.latvec.e31 * ucell.lat0_angstrom;
    dp_cell[7] = ucell.latvec.e32 * ucell.lat0_angstrom;
    dp_cell[8] = ucell.latvec.e33 * ucell.lat0_angstrom;

    dp_coord.resize(3 * ucell.nat);
    const int nat = ucell.nat;
#pragma omp parallel for schedule(static) if (nat >= 256)
    for (int iat = 0; iat < nat; ++iat)
    {
        const int it = atom_type_index[iat];
        const int ia = atom_local_index[iat];
        dp_coord[3 * iat]     = ucell.atoms[it].tau[ia].x * ucell.lat0_angstrom;
        dp_coord[3 * iat + 1] = ucell.atoms[it].tau[ia].y * ucell.lat0_angstrom;
        dp_coord[3 * iat + 2] = ucell.atoms[it].tau[ia].z * ucell.lat0_angstrom;
    }

#ifdef __DPMD
    dp_potential = 0;
    dp_force.zero_out();
    dp_virial.zero_out();
    dp_model_force.clear();
    dp_model_virial.clear();

    dp.compute(dp_potential, dp_model_force, dp_model_virial, dp_coord, atype, dp_cell, fparam, aparam);

    // rescale the energy, force, and stress
    const double fact_e = rescaling / ModuleBase::Ry_to_eV;
    const double fact_f = rescaling / (ModuleBase::Ry_to_eV * ModuleBase::ANGSTROM_AU);
    const double fact_v = rescaling / (ucell.omega * ModuleBase::Ry_to_eV);

    dp_potential *= fact_e;
    GlobalV::ofs_running << " #TOTAL ENERGY# " << std::setprecision(11) << dp_potential * ModuleBase::Ry_to_eV << " eV"
                         << std::endl;

    const int nat_f = ucell.nat;
#pragma omp parallel for schedule(static) if (nat_f >= 256)
    for (int i = 0; i < nat_f; ++i)
    {
        dp_force(i, 0) = dp_model_force[3 * i] * fact_f;
        dp_force(i, 1) = dp_model_force[3 * i + 1] * fact_f;
        dp_force(i, 2) = dp_model_force[3 * i + 2] * fact_f;
    }

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            dp_virial(i, j) = dp_model_virial[3 * i + j] * fact_v;
        }
    }
#else
    ModuleBase::WARNING_QUIT("ESolver_DP", "Please recompile with -D__DPMD");
#endif
    ModuleBase::timer::end("ESolver_DP", "runner");
}

double ESolver_DP::cal_energy()
{
    return dp_potential;
}

void ESolver_DP::cal_force(BaseCell& basecell, ModuleBase::matrix& force)
{
    basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    force = dp_force;
    ModuleIO::print_force(GlobalV::ofs_running, ucell, "TOTAL-FORCE (eV/Angstrom)", force, false);
}

void ESolver_DP::cal_stress(BaseCell& basecell, ModuleBase::matrix& stress)
{
    basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    stress = dp_virial;

    ModuleIO::print_stress("TOTAL-STRESS", stress, true, false, GlobalV::ofs_running);

    // external stress
    double unit_transform = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
    double external_stress[3] = {this->inp_->press1, this->inp_->press2, this->inp_->press3};
    for (int i = 0; i < 3; i++)
    {
        stress(i, i) -= external_stress[i] / unit_transform;
    }
}

void ESolver_DP::after_all_runners(BaseCell& basecell)
{
    basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    GlobalV::ofs_running << "\n --------------------------------------------" << std::endl;
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
            iat++;
        }
    }
    assert(ucell.nat == iat);
}
#endif
