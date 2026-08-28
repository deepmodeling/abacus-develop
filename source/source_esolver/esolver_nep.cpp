/**
 * @file esolver_nep.cpp
#include "source_io/module_parameter/parameter.h"
 * @brief Implementation of ESolver_NEP class for neuroevolution potential (NEP).
 *
 * This file contains the implementation of the ESolver_NEP class, which is used for solving the energy and forces in a
 * NEP simulation.
 * NEP is a method for training deep neural networks to accurately predict the potential energy surface of a
 * molecular system.
 *
 * For more information about NEP, see the following reference:
 * 1. https://gpumd.org/potentials/nep.html
 * 2. https://doi.org/10.1002/mgea.70028
 *
 * @author MoseyQAQ
 * @date 2025-10-10
 */
#include "esolver_nep.h"
#include "source_base/parallel_common.h"
#include "source_base/timer.h"
#include "source_io/module_output/output_log.h"
#include "source_io/module_parameter/parameter.h"

#include <algorithm>
#include <unordered_map>

using namespace ModuleESolver;

void ESolver_NEP::before_all_runners(BaseCell& basecell, const Input_para& inp)
{
    basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    this->inp_ = &inp;

    nep_potential = 0.0;
    nep_force.create(ucell.nat, 3);
    nep_virial.create(3, 3);
    atype.resize(ucell.nat);
    nep_cell.resize(9);
    nep_coord.resize(3 * ucell.nat);
    nep_virial_sum.resize(9);
    _e.resize(ucell.nat);
    _f.resize(3 * ucell.nat);
    _v.resize(9 * ucell.nat);

#ifdef __NEP
    /// determine the type map from STRU to NEP model
    type_map(ucell);
#endif
}

void ESolver_NEP::runner(BaseCell& basecell, const int istep)
{
    basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    ModuleBase::TITLE("ESolver_NEP", "runner");
    ModuleBase::timer::start("ESolver_NEP", "runner");

    // note that NEP are column major, thus a transpose is needed
    // cell
    nep_cell[0] = ucell.latvec.e11 * ucell.lat0_angstrom;
    nep_cell[1] = ucell.latvec.e21 * ucell.lat0_angstrom;
    nep_cell[2] = ucell.latvec.e31 * ucell.lat0_angstrom;
    nep_cell[3] = ucell.latvec.e12 * ucell.lat0_angstrom;
    nep_cell[4] = ucell.latvec.e22 * ucell.lat0_angstrom;
    nep_cell[5] = ucell.latvec.e32 * ucell.lat0_angstrom;
    nep_cell[6] = ucell.latvec.e13 * ucell.lat0_angstrom;
    nep_cell[7] = ucell.latvec.e23 * ucell.lat0_angstrom;
    nep_cell[8] = ucell.latvec.e33 * ucell.lat0_angstrom;

    // coord
    nep_coord.resize(3 * ucell.nat);
    const int nat = ucell.nat;
#pragma omp parallel for schedule(static) if (nat >= 256)
    for (int iat = 0; iat < nat; ++iat)
    {
        const int it = atom_type_index[iat];
        const int ia = atom_local_index[iat];
        nep_coord[iat] = ucell.atoms[it].tau[ia].x * ucell.lat0_angstrom;
        nep_coord[iat + nat] = ucell.atoms[it].tau[ia].y * ucell.lat0_angstrom;
        nep_coord[iat + 2 * nat] = ucell.atoms[it].tau[ia].z * ucell.lat0_angstrom;
    }

#ifdef __NEP
    nep_potential = 0.0;
    nep_force.zero_out();
    nep_virial.zero_out();

    nep.compute(atype, nep_cell, nep_coord, _e, _f, _v);

    // unit conversion
    const double fact_e = 1.0 / ModuleBase::Ry_to_eV;
    const double fact_f = 1.0 / (ModuleBase::Ry_to_eV * ModuleBase::ANGSTROM_AU);
    const double fact_v = 1.0 / (ucell.omega * ModuleBase::Ry_to_eV);

    // potential energy
    double energy_sum = 0.0;
#pragma omp parallel for reduction(+:energy_sum) schedule(static) if (nat >= 256)
    for (int i = 0; i < nat; ++i)
    {
        energy_sum += _e[i];
    }
    nep_potential = fact_e * energy_sum;
    GlobalV::ofs_running << " #TOTAL ENERGY# " << std::setprecision(11) << nep_potential * ModuleBase::Ry_to_eV << " eV"
                         << std::endl;

    // forces
#pragma omp parallel for schedule(static) if (nat >= 256)
    for (int i = 0; i < nat; ++i)
    {
        nep_force(i, 0) = _f[i] * fact_f;
        nep_force(i, 1) = _f[i + nat] * fact_f;
        nep_force(i, 2) = _f[i + 2 * nat] * fact_f;
    }

    // virial
    double v0 = 0.0;
    double v1 = 0.0;
    double v2 = 0.0;
    double v3 = 0.0;
    double v4 = 0.0;
    double v5 = 0.0;
    double v6 = 0.0;
    double v7 = 0.0;
    double v8 = 0.0;
#pragma omp parallel for reduction(+:v0, v1, v2, v3, v4, v5, v6, v7, v8) schedule(static) if (nat >= 256)
    for (int i = 0; i < nat; ++i)
    {
        v0 += _v[i];
        v1 += _v[nat + i];
        v2 += _v[2 * nat + i];
        v3 += _v[3 * nat + i];
        v4 += _v[4 * nat + i];
        v5 += _v[5 * nat + i];
        v6 += _v[6 * nat + i];
        v7 += _v[7 * nat + i];
        v8 += _v[8 * nat + i];
    }
    nep_virial_sum[0] = v0;
    nep_virial_sum[1] = v1;
    nep_virial_sum[2] = v2;
    nep_virial_sum[3] = v3;
    nep_virial_sum[4] = v4;
    nep_virial_sum[5] = v5;
    nep_virial_sum[6] = v6;
    nep_virial_sum[7] = v7;
    nep_virial_sum[8] = v8;

    // virial -> stress
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            nep_virial(i, j) = nep_virial_sum[3 * i + j] * fact_v;
        }
    }
#else
    ModuleBase::WARNING_QUIT("ESolver_NEP", "Please recompile with -D__NEP");
#endif
    ModuleBase::timer::end("ESolver_NEP", "runner");
}

double ESolver_NEP::cal_energy()
{
    return nep_potential;
}

void ESolver_NEP::cal_force(BaseCell& basecell, ModuleBase::matrix& force)
{
    basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    force = nep_force;
    ModuleIO::print_force(GlobalV::ofs_running, ucell, "TOTAL-FORCE (eV/Angstrom)", force, false);
}

void ESolver_NEP::cal_stress(BaseCell& basecell, ModuleBase::matrix& stress)
{
    basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    stress = nep_virial;
    ModuleIO::print_stress("TOTAL-STRESS", stress, true, false, GlobalV::ofs_running);

    // external stress
    double unit_transform = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
    double external_stress[3] = {this->inp_->press1, this->inp_->press2, this->inp_->press3};
    for (int i = 0; i < 3; i++)
    {
        stress(i, i) -= external_stress[i] / unit_transform;
    }
}

void ESolver_NEP::after_all_runners(BaseCell& basecell)
{
    basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    GlobalV::ofs_running << "\n --------------------------------------------" << std::endl;
    GlobalV::ofs_running << std::setprecision(16);
    GlobalV::ofs_running << " !FINAL_ETOT_IS " << nep_potential * ModuleBase::Ry_to_eV << " eV" << std::endl;
    GlobalV::ofs_running << " --------------------------------------------\n\n" << std::endl;
}

#ifdef __NEP
void ESolver_NEP::type_map(const UnitCell& ucell)
{
    // parse the element list from NEP model file
    std::unordered_map<std::string, int> label;
    std::string temp;
    for (int i = 0; i < nep.element_list.size(); ++i)
    {
        label[nep.element_list[i]] = i; //> label: map from element string to index int.
    }

    std::cout << "\n Element list of model file " << nep_file << " " << std::endl;
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

    // parse the atype based on the element list
    int iat = 0;
    for (int it = 0; it < ucell.ntype; ++it)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
        {
            if (label.find(ucell.atoms[it].label) == label.end())
            {
                ModuleBase::WARNING_QUIT("ESolver_NEP",
                                         "The label " + ucell.atoms[it].label + " is not found in the type map.");
            }
            atype[iat] = label[ucell.atoms[it].label];
            iat++;
        }
    }
    assert(ucell.nat == iat);
}
#endif
