#include "check_atomic_stru.h"

#include "module_base/element_covalent_radius.h"

void Check_Atomic_Stru::check_atomic_stru(UnitCell& ucell, const double& factor) {
    // First we calculate all bond length in the structure,
    // and compare with the covalent_bond_length,
    // if there has bond length is shorter than covalent_bond_length * factor,
    // we think this structure is unreasonable.
    const double warning_coef = 0.6;
    assert(ucell.ntype > 0);
    std::stringstream errorlog;
    bool all_pass = true;
    bool no_warning = true;
    for (int it1 = 0; it1 < ucell.ntype; it1++) {
        std::string symbol1 = "";
        for (char ch: ucell.atoms[it1].label) {
            if (std::isalpha(ch)) {
                symbol1.push_back(ch);
            }
        }
        // std::string symbol1 = ucell.atoms[it1].label;
        double symbol1_covalent_radius;
        if (ModuleBase::CovalentRadius.find(symbol1)
            != ModuleBase::CovalentRadius.end()) {
            symbol1_covalent_radius = ModuleBase::CovalentRadius.at(symbol1);
        } else {
            std::stringstream mess;
            mess << "Notice: symbol '" << symbol1
                 << "' is not an element symbol!!!! ";
            mess << "set the covalent radius to be 0." << std::endl;
            GlobalV::ofs_running << mess.str();
            std::cout << mess.str();
            symbol1_covalent_radius = 0.0;
        }

        for (int ia1 = 0; ia1 < ucell.atoms[it1].na; ia1++) {
            double x1 = ucell.atoms[it1].taud[ia1].x;
            double y1 = ucell.atoms[it1].taud[ia1].y;
            double z1 = ucell.atoms[it1].taud[ia1].z;

            for (int it2 = 0; it2 < ucell.ntype; it2++) {
                std::string symbol2 = ucell.atoms[it2].label;
                double symbol2_covalent_radius;
                if (ModuleBase::CovalentRadius.find(symbol2)
                    != ModuleBase::CovalentRadius.end()) {
                    symbol2_covalent_radius
                        = ModuleBase::CovalentRadius.at(symbol2);
                } else {
                    symbol2_covalent_radius = 0.0;
                }

                double covalent_length
                    = (symbol1_covalent_radius + symbol2_covalent_radius)
                      / ModuleBase::BOHR_TO_A;

                for (int ia2 = 0; ia2 < ucell.atoms[it2].na; ia2++) {
                    for (int a = -1; a < 2; a++) {
                        for (int b = -1; b < 2; b++) {
                            for (int c = -1; c < 2; c++) {
                                if (it1 > it2) {
                                    continue;
                                } else if (it1 == it2 && ia1 > ia2) {
                                    continue;
                                } else if (it1 == it2 && ia1 == ia2 && a == 0
                                         && b == 0 && c == 0) {
                                    continue;
}

                                double x2 = ucell.atoms[it2].taud[ia2].x + a;
                                double y2 = ucell.atoms[it2].taud[ia2].y + b;
                                double z2 = ucell.atoms[it2].taud[ia2].z + c;

                                double bond_length
                                    = sqrt(pow((x2 - x1) * ucell.a1.x
                                                   + (y2 - y1) * ucell.a2.x
                                                   + (z2 - z1) * ucell.a3.x,
                                               2)
                                           + pow((x2 - x1) * ucell.a1.y
                                                     + (y2 - y1) * ucell.a2.y
                                                     + (z2 - z1) * ucell.a3.y,
                                                 2)
                                           + pow((x2 - x1) * ucell.a1.z
                                                     + (y2 - y1) * ucell.a2.z
                                                     + (z2 - z1) * ucell.a3.z,
                                                 2))
                                      * ucell.lat0;

                                if (bond_length < covalent_length * factor
                                    || bond_length
                                           < covalent_length * warning_coef) {
                                    errorlog.setf(std::ios_base::fixed,
                                                  std::ios_base::floatfield);
                                    errorlog << std::setw(3) << ia1 + 1
                                             << "-th " << std::setw(3)
                                             << ucell.atoms[it1].label << ", ";
                                    errorlog << std::setw(3) << ia2 + 1
                                             << "-th " << std::setw(3)
                                             << ucell.atoms[it2].label;
                                    errorlog << " (cell:" << std::setw(2) << a
                                             << " " << std::setw(2) << b << " "
                                             << std::setw(2) << c << ")";
                                    errorlog << ", distance= "
                                             << std::setprecision(3)
                                             << bond_length << " Bohr (";
                                    errorlog
                                        << bond_length * ModuleBase::BOHR_TO_A
                                        << " Angstrom)" << std::endl;

                                    if (bond_length
                                        < covalent_length * factor) {
                                        all_pass = false;
                                    } else {
                                        no_warning = false;
                                    }
                                }
                            } // c
                        }     // b
                    }         // a
                }             // ia2
            }                 // it2
        }                     // ia1
    }                         // it1

    if (!all_pass || !no_warning) {
        std::stringstream mess;
        mess << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
             << std::endl;
        mess << "%%%%%% WARNING  WARNING  WARNING  WARNING  WARNING  %%%%%%"
             << std::endl;
        mess << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
             << std::endl;
        mess << "!!! WARNING: Some atoms are too close!!!" << std::endl;
        mess << "!!! Please check the nearest-neighbor list in log file."
             << std::endl;
        mess << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
             << std::endl;
        mess << "%%%%%% WARNING  WARNING  WARNING  WARNING  WARNING  %%%%%%"
             << std::endl;
        mess << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
             << std::endl;

        GlobalV::ofs_running << mess.str() << mess.str() << mess.str()
                             << errorlog.str();
        std::cout << mess.str() << mess.str() << mess.str() << std::endl;

        if (!all_pass) {
            mess.clear();
            mess.str("");
            mess << "If this structure is what you want, you can set "
                    "'min_dist_coef'"
                 << std::endl;
            mess << "as a smaller value (the current value is " << factor
                 << ") in INPUT file." << std::endl;
            GlobalV::ofs_running << mess.str();
            std::cout << mess.str();
            ModuleBase::WARNING_QUIT("Input", "The structure is unreasonable!");
        }
    }
}
