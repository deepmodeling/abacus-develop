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

namespace ModuleESolver
{

    void ESolver_DP::init(Input& inp, UnitCell& ucell)
    {
        ucell_ = &ucell;
        dp_potential = 0;
        dp_force.create(ucell.nat, 3);
        dp_virial.create(3, 3);

        cell.resize(9);
        atype.resize(ucell.nat);
        coord.resize(3 * ucell.nat);
        dp_type.resize(ucell.ntype);

        bool find_type = type_map(ucell);

        /// determine the type map from STRU to DP model
        int iat = 0;
        for (int it = 0; it < ucell.ntype; ++it)
        {
            for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
            {
                if (find_type)
                {
                    atype[iat] = dp_type[it];
                }
                else
                {
                    atype[iat] = it;
                }
                iat++;
            }
        }
        assert(ucell.nat == iat);
    }

    void ESolver_DP::run(const int istep, UnitCell& ucell)
    {
        ModuleBase::TITLE("ESolver_DP", "Run");
        ModuleBase::timer::tick("ESolver_DP", "Run");

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

        dp.compute(dp_potential, f, v, coord, atype, cell);

        dp_potential /= ModuleBase::Ry_to_eV;
        GlobalV::ofs_running << " final etot is " << std::setprecision(11) << dp_potential * ModuleBase::Ry_to_eV
                             << " eV" << std::endl;

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
        ModuleBase::timer::tick("ESolver_DP", "Run");
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
        double external_stress[3] = {GlobalV::PRESS1, GlobalV::PRESS2, GlobalV::PRESS3};
        for (int i = 0; i < 3; i++)
        {
            stress(i, i) -= external_stress[i] / unit_transform;
        }

        ModuleIO::print_stress("TOTAL-STRESS", stress, true, false);
    }

    void ESolver_DP::post_process(void)
    {
        GlobalV::ofs_running << "\n\n --------------------------------------------" << std::endl;
        GlobalV::ofs_running << std::setprecision(16);
        GlobalV::ofs_running << " !FINAL_ETOT_IS " << dp_potential * ModuleBase::Ry_to_eV << " eV" << std::endl;
        GlobalV::ofs_running << " --------------------------------------------\n\n" << std::endl;
    }

    bool ESolver_DP::type_map(const UnitCell& ucell)
    {
        bool ok = false;
        bool find_type = false;

        if (GlobalV::MY_RANK == 0)
        {
            std::ifstream ifs(dp_file);
            std::string word;
            int ntype_dp = 0;
            std::string* label = nullptr;

            if (ifs)
            {
                ok = true;
            }

            if (ok)
            {
                while (std::getline(ifs, word, '"'))
                {
                    if (word == "type_map")
                    {
                        find_type = true;
                        break;
                    }
                }

                if (find_type)
                {
                    std::getline(ifs, word, '"'); ///< the string is ":[", which is useless
                    std::stringstream ss;
                    while (word[0] != ']')
                    {
                        std::getline(ifs, word, '"');
                        ++ntype_dp;
                        ss << word << "  ";
                        std::getline(ifs, word, '"');
                    }

                    GlobalV::ofs_running << std::endl;
                    GlobalV::ofs_running << "Determine the type map from DP model" << std::endl;
                    GlobalV::ofs_running << "ntype read from DP model: " << ntype_dp << std::endl;

                    label = new std::string[ntype_dp];
                    for (int it = 0; it < ntype_dp; ++it)
                    {
                        ss >> label[it];
                        GlobalV::ofs_running << "  " << label[it];
                    }
                    GlobalV::ofs_running << std::endl << std::endl;

                    for (int it = 0; it < ucell.ntype; ++it)
                    {
                        bool consistent = false;
                        for (int it2 = 0; it2 < ntype_dp; ++it2)
                        {
                            if (ucell.atom_label[it] == label[it2])
                            {
                                dp_type[it] = it2;
                                consistent = true;
                            }
                        }
                        if (!consistent)
                        {
                            ModuleBase::WARNING_QUIT("ESolver_DP", "Unsupported atom types for the DP model");
                        }
                    }
                    delete[] label;
                }
            }
            ifs.close();
        }

#ifdef __MPI
        Parallel_Common::bcast_bool(ok);
        Parallel_Common::bcast_bool(find_type);
        Parallel_Common::bcast_int(dp_type.data(), ucell.ntype);
#endif

        if (!ok)
        {
            ModuleBase::WARNING_QUIT("ESolver_DP", "can not find the DP model");
        }
        return find_type;
    }
}
