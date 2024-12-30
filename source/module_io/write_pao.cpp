#include "write_pao.h"
#include "module_base/global_variable.h"
#include "module_parameter/parameter.h"

#include <fstream>
namespace ModuleIO
{
void print_PAOs(const UnitCell& ucell)
{
    if (GlobalV::MY_RANK != 0)
    {
        return;
    }
    for (int it = 0; it < ucell.ntype; it++)
    {
        for (int icc = 0; icc < ucell.atoms[it].ncpp.nchi; icc++)
        {
            std::stringstream ss;
            ss << PARAM.globalv.global_out_dir << ucell.atoms[it].label << "/" << ucell.atoms[it].label << "-"
               << ucell.atoms[it].ncpp.els[icc] << ".ORBITAL";

            std::ofstream ofs(ss.str().c_str());
            ofs << "Mesh " << ucell.atoms[it].ncpp.msh;
            ofs << "\n" << std::setw(15) << "Radial" << std::setw(15) << "Psi" << std::setw(15) << "Rab";

            for (int i = 0; i < ucell.atoms[it].ncpp.msh; i++)
            {
                ofs << "\n"
                    << std::setw(15) << ucell.atoms[it].ncpp.r[i] << std::setw(15) << ucell.atoms[it].ncpp.chi(icc, i)
                    << std::setw(15) << ucell.atoms[it].ncpp.rab[i];
            }
            ofs.close();
        }
        // end out put
    }
    return;
}
} // namespace ModuleIO