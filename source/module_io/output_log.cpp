#include "output_log.h"

#include "module_base/constants.h"
#include "module_base/formatter.h"
#include "module_base/global_variable.h"

namespace ModuleIO
{
void output_convergence_after_scf(bool& convergence, double& energy, std::ofstream& ofs_running)
{
    if (convergence)
    {
        ofs_running << "\n charge density convergence is achieved" << std::endl;
        ofs_running << " final etot is " << std::setprecision(11) << energy * ModuleBase::Ry_to_eV << " eV" << std::endl;
    }
    else
    {
        ofs_running << " !! convergence has not been achieved @_@" << std::endl;
        std::cout << " !! CONVERGENCE HAS NOT BEEN ACHIEVED !!" << std::endl;
    }
}

void output_efermi(bool& convergence, double& efermi, std::ofstream& ofs_running)
{
    if (convergence && GlobalV::OUT_LEVEL != "m")
    {
        ofs_running << std::setprecision(16);
        ofs_running << " EFERMI = " << std::setprecision(11) << efermi * ModuleBase::Ry_to_eV << " eV" << std::endl;
    }
}

void print_force(std::ofstream& ofs_running,
                 const UnitCell& cell,
                 const std::string& name,
                 const ModuleBase::matrix& force,
                 bool ry)
{
    const double output_acc = 1.0e-8;
    double fac = 1.0;
    if (!ry)
    {
        fac = ModuleBase::Ry_to_eV / 0.529177;
    }

    std::vector<std::string> atom_label;
    std::vector<double> force_x;
    std::vector<double> force_y;
    std::vector<double> force_z;
    std::string table;
    context.set_context({"title", "force", "force", "force"});
    int iat = 0;
    for (int it = 0; it < cell.ntype; it++)
    {
        for (int ia = 0; ia < cell.atoms[it].na; ia++)
        {
            std::string atom_labels = cell.atoms[it].label + std::to_string(ia + 1);
            double fx = std::abs(force(iat, 0)) > output_acc ? force(iat, 0) * fac : 0.0;
            double fy = std::abs(force(iat, 1)) > output_acc ? force(iat, 1) * fac : 0.0;
            double fz = std::abs(force(iat, 2)) > output_acc ? force(iat, 2) * fac : 0.0;
            atom_label.push_back(atom_labels);
            force_x.push_back(fx);
            force_y.push_back(fy);
            force_z.push_back(fz);

            iat++;
        }
    }

    context.enable_title();
    context << "atom" << atom_label << "x" << force_x << "y" << force_y << "z" << force_z;
    context.center_title();
    context.set_overall_title("TOTAL-FORCE (eV/Angstrom)");
    table = context.str();
    ofs_running << table << std::endl;
    if (GlobalV::TEST_FORCE)
    {
        std::cout << table << std::endl;
    }
    return;
}

void print_stress(const std::string& name, const ModuleBase::matrix& scs, const bool screen, const bool ry)
{
    const double output_acc = 1.0e-8;
    double unit_transform = 1;
    std::string title = name;
    if (ry)
    {
        title += " (a.u.)";
    }
    else
    {
        title += " (KBAR)";
        unit_transform = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
    }

    std::vector<double> stress_x;
    std::vector<double> stress_y;
    std::vector<double> stress_z;
    std::string table;
    context.set_context({"double_w20_f10", "double_w20_f10", "double_w20_f10"});
    for (int i = 0; i < 3; i++)
    {
        double sx = scs(i, 0) * unit_transform;
        double sy = scs(i, 1) * unit_transform;
        double sz = scs(i, 2) * unit_transform;
        stress_x.push_back(sx);
        stress_y.push_back(sy);
        stress_z.push_back(sz);
    }

    context.enable_title();
    context << title.c_str() << stress_x << " " << stress_y << " " << stress_z;
    context.center_title();
    table = context.str();
    GlobalV::ofs_running << table << std::endl;
    if (screen)
    {
        std::cout << table << std::endl;
    }
    return;
}

}// namespace ModuleIO