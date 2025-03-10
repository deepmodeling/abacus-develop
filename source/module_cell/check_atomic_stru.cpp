#include "check_atomic_stru.h"

#include "module_base/element_covalent_radius.h"
#include "module_base/timer.h"

namespace unitcell
{

void check_atomic_stru(UnitCell& ucell, const double& factor)
{
    ModuleBase::timer::tick("unitcell", "check_atomic_stru");
    // First we calculate all bond length in the structure,
    // and compare with the covalent_bond_length,
    // if there has bond length is shorter than covalent_bond_length * factor,
    // we think this structure is unreasonable.
    assert(ucell.ntype > 0);
    bool all_pass = true;
    bool no_warning = true;
    std::stringstream errorlog;
    errorlog.setf(std::ios_base::fixed, std::ios_base::floatfield);

    if (GlobalV::MY_RANK == 0)
    {
        
        const int ntype = ucell.ntype;
        const double lat0 = ucell.lat0;
        const double warning_coef = 0.6;
        const double max_factor_coef = std::max(warning_coef, factor);

        std::vector<double> symbol_covalent_radiuss(ntype);
        for (int it = 0; it < ntype; it++)
        {
            std::string symbol1 = "";
            for (char ch: ucell.atoms[it].label)
            {
                if (std::isalpha(ch))
                {
                    symbol1.push_back(ch);
                }
            }

            if (ModuleBase::CovalentRadius.find(symbol1) != ModuleBase::CovalentRadius.end())
            {
                symbol_covalent_radiuss[it] = ModuleBase::CovalentRadius.at(symbol1);
            }
            else
            {
                std::stringstream mess;
                mess << "Notice: symbol '" << symbol1 << "' is not an element symbol!!!! ";
                mess << "set the covalent radius to be 0." << std::endl;
                GlobalV::ofs_running << mess.str();
                std::cout << mess.str();
            }
        }
        std::vector<double> latvec (9);
        latvec[0] = ucell.a1.x;
        latvec[1] = ucell.a2.x;
        latvec[2] = ucell.a3.x;
        latvec[3] = ucell.a1.y;
        latvec[4] = ucell.a2.y;
        latvec[5] = ucell.a3.y;
        latvec[6] = ucell.a1.z;
        latvec[7] = ucell.a2.z;
        latvec[8] = ucell.a3.z;
        std::vector<double> A(27*3);
        std::vector<std::string> cell(27);
        std::vector<std::string> label(ntype);
        for (int i = 0; i < 27; i++)
        {
            int a = (i / 9) % 3 - 1;
            int b = (i / 3) % 3 - 1;
            int c = i % 3 - 1;
            A[3 * i] = a * latvec[0] + b * latvec[1] + c * latvec[2];
            A[3 * i + 1] = a * latvec[3] + b * latvec[4] + c * latvec[5];
            A[3 * i + 2] = a * latvec[6] + b * latvec[7] + c * latvec[8];
            std::ostringstream tmp_oss;
            tmp_oss << " (cell:" << std::setw(2) << a << " " << std::setw(2) << b << " " << std::setw(2) << c
                    << "), distance= ";
            cell[i] = tmp_oss.str();
        }
        for (int it = 0; it < ntype; it++)
        {
            std::ostringstream tmp_oss;
            tmp_oss << std::setw(3) << ucell.atoms[it].label;
            label[it] = tmp_oss.str();
        }

        const double bohr_to_a = ModuleBase::BOHR_TO_A;
#pragma omp parallel
        {
            std::vector<double> delta_lat(3);
#pragma omp for schedule(dynamic)
            for (int iat = 0; iat < ucell.nat; iat++)
            {
                const int it1 = ucell.iat2it[iat];
                const int ia1 = ucell.iat2ia[iat];
                const double symbol1_covalent_radius = symbol_covalent_radiuss[it1];
                double x1 = ucell.atoms[it1].taud[ia1].x;
                double y1 = ucell.atoms[it1].taud[ia1].y;
                double z1 = ucell.atoms[it1].taud[ia1].z;
                for (int it2 = it1; it2 < ntype; it2++)
                {
                    double symbol2_covalent_radius = symbol_covalent_radiuss[it2];
                    double covalent_length = (symbol1_covalent_radius + symbol2_covalent_radius) / bohr_to_a;
                    const double max_error = covalent_length * max_factor_coef / ucell.lat0;
                    const double max_error_2 = max_error * max_error;
                    const double factor_error = covalent_length * factor;
                    for (int ia2 = ia1; ia2 < ucell.atoms[it2].na; ia2++)
                    {
                        const bool is_same_atom = (it1 == it2) && (ia1 == ia2);
                        double delta_x = ucell.atoms[it2].taud[ia2].x - x1;
                        double delta_y = ucell.atoms[it2].taud[ia2].y - y1;
                        double delta_z = ucell.atoms[it2].taud[ia2].z - z1;
                        delta_lat[0] = delta_x * latvec[0] + delta_y * latvec[1] + delta_z * latvec[2];
                        delta_lat[1] = delta_x * latvec[3] + delta_y * latvec[4] + delta_z * latvec[5];
                        delta_lat[2] = delta_x * latvec[6] + delta_y * latvec[7] + delta_z * latvec[8];
                        for (int i = 0; i < 27; i++)
                        {
                            if ((is_same_atom) && (i == 13))
                                continue;
                            const int offset = i * 3;
                            const double part1 = delta_lat[0] + A[offset];
                            const double part2 = delta_lat[1] + A[offset + 1];
                            const double part3 = delta_lat[2] + A[offset + 2];
                            const double bond_length = part1 * part1 + part2 * part2 + part3 * part3;
                            const bool flag = bond_length < max_error_2 ? true : false;
                            if (flag)
                            {
                                const double sqrt_bon = sqrt(bond_length) * lat0;
                                #pragma omp critical
                                {
                                    no_warning = false;
                                    all_pass = all_pass && (sqrt_bon < factor_error ? false : true);
                                    errorlog << std::setw(3) << ia1 + 1 << "-th " << label[it1] << ", " << std::setw(3)
                                             << ia2 + 1 << "-th " << label[it2] << cell[i] << std::setprecision(3)
                                             << sqrt_bon << " Bohr (" << sqrt_bon * bohr_to_a << " Angstrom)\n";
                                }
                            }
                        } 
                    } // ia2
                } // it2
            } // iat
        }
    }
    if (!all_pass || !no_warning)
    {
        std::stringstream mess;
        mess << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
        mess << "%%%%%% WARNING  WARNING  WARNING  WARNING  WARNING  %%%%%%" << std::endl;
        mess << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
        mess << "!!! WARNING: Some atoms are too close!!!" << std::endl;
        mess << "!!! Please check the nearest-neighbor list in log file." << std::endl;
        mess << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
        mess << "%%%%%% WARNING  WARNING  WARNING  WARNING  WARNING  %%%%%%" << std::endl;
        mess << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;

        GlobalV::ofs_running << mess.str() << mess.str() << mess.str() << errorlog.str();
        std::cout << mess.str() << mess.str() << mess.str() << std::endl;
        if (!all_pass)
        {
            mess.clear();
            mess.str("");
            mess << "If this structure is what you want, you can set 'min_dist_coef'\n";
            mess << "as a smaller value (the current value is " << factor << ") in INPUT file." << std::endl;
            GlobalV::ofs_running << mess.str();
            std::cout << mess.str();
            ModuleBase::WARNING_QUIT("Input", "The structure is unreasonable!");
        }
    }

    ModuleBase::timer::tick("unitcell", "check_atomic_stru");
}

}
