#include "check_atomic_stru.h"
#include "module_base/timer.h"
#include "module_base/element_covalent_radius.h"

void Check_Atomic_Stru::check_atomic_stru(UnitCell& ucell, const double& factor) {
    // First we calculate all bond length in the structure,
    // and compare with the covalent_bond_length,
    // if there has bond length is shorter than covalent_bond_length * factor,
    // we think this structure is unreasonable.
    ModuleBase::timer::tick("Check_Atomic_Stru", "Check_Atomic_Stru");
    assert(ucell.ntype > 0);
    const int ntype = ucell.ntype;
    const double lat0 = ucell.lat0;
    bool all_pass = true;
    bool no_warning = true;
    const double warning_coef = 0.6;
    const double min_factor_coef=std::min(warning_coef,factor);

    std::stringstream errorlog;
    errorlog.setf(std::ios_base::fixed, std::ios_base::floatfield);

    std::vector<double> symbol_covalent_radiuss(ntype);
    for (int it = 0;it < ntype ; it++)
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
            symbol_covalent_radiuss[it] =  ModuleBase::CovalentRadius.at(symbol1);
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
    double * latvec = new double [9];
    latvec[0]=ucell.a1.x; latvec[1]=ucell.a2.x; latvec[2]=ucell.a3.x;
    latvec[3]=ucell.a1.y; latvec[4]=ucell.a2.y; latvec[5]=ucell.a3.y;
    latvec[6]=ucell.a1.z; latvec[7]=ucell.a2.z; latvec[8]=ucell.a3.z;
    std::vector<double> A(27*3);
    for (int i=0;i<27;i++)
    {
        int a = (i / 9) % 3 - 1; // 计算a的值
        int b = (i / 3) % 3 - 1; // 计算b的值
        int c = i % 3 - 1;       // 计算c的值
        A[3*i]= a*latvec[0] + b *latvec[1] + c* latvec[2];
        A[3*i+1]= a*latvec[3] + b *latvec[4] + c* latvec[5];
        A[3*i+2]= a*latvec[6] + b *latvec[7] + c* latvec[8];
    }
    for (int it1 = 0; it1 < ntype; it1++) {
        const double symbol1_covalent_radius=symbol_covalent_radiuss[it1];
        for (int ia1 = 0; ia1 < ucell.atoms[it1].na; ia1++) {
            double x1 = ucell.atoms[it1].taud[ia1].x;
            double y1 = ucell.atoms[it1].taud[ia1].y;
            double z1 = ucell.atoms[it1].taud[ia1].z;

            for (int it2 = it1; it2 < ucell.ntype; it2++) {
                std::string symbol2 = ucell.atoms[it2].label;
                double symbol2_covalent_radius=symbol_covalent_radiuss[it1];

                double covalent_length = (symbol1_covalent_radius + symbol2_covalent_radius) / ModuleBase::BOHR_TO_A;
                const double min_error = covalent_length * min_factor_coef / ucell.lat0;
                const double min_error_2 = min_error * min_error;

                for (int ia2 = ia1; ia2 < ucell.atoms[it2].na; ia2++) 
                {
                    double delta_x = ucell.atoms[it2].taud[ia2].x - x1;
                    double delta_y = ucell.atoms[it2].taud[ia2].y - y1;
                    double delta_z = ucell.atoms[it2].taud[ia2].z - z1;
                    const double delta_x_lat = delta_x*latvec[0] + delta_y*latvec[1]+ delta_z*latvec[2];
                    const double delta_y_lat = delta_x*latvec[3] + delta_y*latvec[4]+ delta_z*latvec[5];
                    const double delta_z_lat = delta_x*latvec[6] + delta_y*latvec[7]+ delta_z*latvec[8];
                    for (int i=0;i<27;i++)
                    {

                        if (it1 == it2 && ia1 == ia2 && i==13) 
                        {
                            continue;
                        }
                        double part1= delta_x_lat+A[3*i];
                        double part2= delta_y_lat+A[3*i+1];
                        double part3= delta_z_lat+A[3*i+2];
                        double bond_length = part1*part1 + part2*part2+ part3*part3;

                        if (bond_length < min_error_2) 
                        {
                            errorlog << std::setw(3) << ia1 + 1
                                        << "-th " << std::setw(3)
                                        << ucell.atoms[it1].label << ", "
                                        << std::setw(3) << ia2 + 1
                                        << "-th " << std::setw(3)
                                        << ucell.atoms[it2].label
                                        << " (cell:" << std::setw(2) << (i / 9) % 3 - 1
                                        << " " << std::setw(2) << (i / 3) % 3 - 1 << " "
                                        << std::setw(2) << i % 3 - 1 << "), distance= "
                                        << std::setprecision(3)
                                        << sqrt(bond_length)/lat0 << " Bohr ("
                                    << sqrt(bond_length) /lat0 * ModuleBase::BOHR_TO_A
                                    << " Angstrom)\n";
                            if (all_pass || no_warning)
                            {
                                if (bond_length < pow(covalent_length * factor/lat0,2)) 
                                {
                                    all_pass = false;
                                } else {
                                    no_warning = false;
                                }
                            }
                            
                        }
                    }         // a
                }             // ia2
            }                 // it2
        }                     // ia1
    }                         // it1
    delete[] latvec;
    ModuleBase::timer::tick("Check_Atomic_Stru", "Check_Atomic_Stru");
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
        std::ofstream outFile("example2.txt");
        if (outFile.is_open()) {

            // 将stringstream中的内容写入到文件
            outFile << errorlog.str();

            // 关闭文件
            outFile.close();
            std::cout << "数据成功写入到文件。\n";
        } else {
            std::cerr << "无法打开文件进行写入。\n";
        }
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
            // ModuleBase::WARNING_QUIT("Input", "The structure is unreasonable!");
        }
    }

}
