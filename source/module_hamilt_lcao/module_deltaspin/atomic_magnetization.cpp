#include "atomic_magnetization.h"

#include <iostream>

#include "module_base/matrix.h"
#include "module_base/name_angular.h"
#include "module_base/tool_title.h"
#include "module_io/mulliken_charge.h"

void calculate_MW_from_lambda(const int& step,
                              LCAO_Hamilt& uhm,
                              Local_Orbital_Charge& loc,
                              const K_Vectors& kv,
                              const UnitCell& ucell)
{
    ModuleBase::TITLE("module_deltaspin", "calculate_MW_from_lambda");

    ModuleBase::matrix orbMulP;
    orbMulP = ModuleIO::cal_mulliken_k(loc.dm_k, uhm, kv);
    
    std::vector<std::vector<std::vector<double>>> AorbMulP = ModuleIO::convert(orbMulP);
    
    if (GlobalV::MY_RANK == 0)
    {
        const int nlocal = (GlobalV::NSPIN == 4) ? GlobalV::NLOCAL / 2 : GlobalV::NLOCAL;
        std::stringstream as;
        as << GlobalV::global_out_dir << "atomic_magnetization.txt";
        std::ofstream os;
        if (step == 0)
        {
            os.open(as.str().c_str());
        }
        else
        {
            os.open(as.str().c_str(), std::ios::app);
        }
        os << "STEP: " << step << std::endl;
        os << "CALCULATE THE MULLIkEN ANALYSIS OF MAGNETIZATION FOR EACH ATOM" << std::endl;

        double sch = 0.0;
        os << std::setprecision(4);
        for (size_t is = 0; is != GlobalV::NSPIN; ++is)
        {
            if (GlobalV::NSPIN == 4 && is > 0)
                continue;
            double sss = 0.0;
            for (size_t iw = 0; iw != nlocal; ++iw)
            {
                sch += orbMulP(is, iw);
                sss += orbMulP(is, iw);
            }
        }
        os << " Total charge:\t" << sch << std::endl;
        os << "Decomposed Mulliken populations" << std::endl;

        for (size_t i = 0; i != ucell.nat; ++i)
        {
            double total_charge = 0.0, atom_mag = 0.0;
            std::vector<double> total_charge_soc(GlobalV::NSPIN);
            const int t = ucell.iat2it[i];
            int num = 0;
            os << i << std::setw(25) << "Zeta of " << ucell.atoms[t].label << std::setw(30) << "Spin 1" << std::setw(30)
               << "Spin 2" << std::setw(30) << "Spin 3" << std::setw(30) << "Spin 4" << std::endl;

            for (size_t L = 0; L != ucell.atoms[t].nwl + 1; ++L)
            {
                std::vector<double> sum_l(GlobalV::NSPIN, 0.0);
                for (size_t Z = 0; Z != ucell.atoms[t].l_nchi[L]; ++Z)
                {
                    std::vector<double> sum_m(GlobalV::NSPIN, 0.0);
                    for (size_t M = 0; M != (2 * L + 1); ++M)
                    {
                        double spin[4];
                        for (int j = 0; j < 4; j++)
                        {
                            spin[j] = ModuleIO::output_cut(AorbMulP[j][i][num]);
                            sum_m[j] += AorbMulP[j][i][num];
                        }
                        os << ModuleBase::Name_Angular[L][M] << std::setw(25) << Z << std::setw(32) << spin[0]
                           << std::setw(30) << spin[1] << std::setw(30) << spin[2] << std::setw(30) << spin[3]
                           << std::endl;
                        num++;
                    }

                    double spin[4];
                    for (int j = 0; j < 4; j++)
                    {
                        spin[j] = ModuleIO::output_cut(sum_m[j]);
                        sum_l[j] += sum_m[j];
                    }
                    os << "  sum over m " << std::setw(45) << spin[0] << std::setw(30) << spin[1] << std::setw(30)
                       << spin[2] << std::setw(30) << spin[3] << std::endl;
                }

                if (ucell.atoms[t].l_nchi[L])
                {
                    double spin[4];
                    for (int j = 0; j < 4; j++)
                    {
                        spin[j] = ModuleIO::output_cut(sum_l[j]);
                        total_charge_soc[j] += sum_l[j];
                    }
                    os << "  sum over m+zeta " << std::setw(40) << spin[0] << std::setw(30) << spin[1] << std::setw(30)
                       << spin[2] << std::setw(30) << spin[3] << std::endl;
                }
            }

            double spin1 = ModuleIO::output_cut(total_charge_soc[0]);
            double spin2 = ModuleIO::output_cut(total_charge_soc[1]);
            double spin3 = ModuleIO::output_cut(total_charge_soc[2]);
            double spin4 = ModuleIO::output_cut(total_charge_soc[3]);
            os << "Total Charge on atom:  " << ucell.atoms[t].label << std::setw(20) << spin1 << std::endl;
            os << "Total Magnetism on atom:  " << ucell.atoms[t].label << std::setw(20) << "(" << spin2 << ", " << spin3
               << ", " << spin4 << ")" << std::endl;
            os << std::endl << std::endl;
        }
        os.close();
    }
}