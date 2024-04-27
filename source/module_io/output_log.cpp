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

    std::string table;
    int iat = 0;
    for (int it = 0; it < cell.ntype; it++)
    {
        for (int ia = 0; ia < cell.atoms[it].na; ia++)
        {
            std::string atom_labels = cell.atoms[it].label + std::to_string(ia + 1);
            double fx = std::abs(force(iat, 0)) > output_acc ? force(iat, 0) * fac : 0.0;
            double fy = std::abs(force(iat, 1)) > output_acc ? force(iat, 1) * fac : 0.0;
            double fz = std::abs(force(iat, 2)) > output_acc ? force(iat, 2) * fac : 0.0;
            table += FmtCore::format("%-10s %20.10f %20.10f %20.10f\n", atom_labels, fx, fy, fz);
            iat++;
        }
    }
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
    std::string unit = "";
    if (ry)
    {
        title += " (a.u.)";
        unit = " a.u.";
    }
    else
    {
        title += " (KBAR)";
        unit = " KBAR";
        unit_transform = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
    }

    std::string table;
    for (int i = 0; i < 3; i++)
    {
        double sx = scs(i, 0) * unit_transform;
        double sy = scs(i, 1) * unit_transform;
        double sz = scs(i, 2) * unit_transform;
        table += FmtCore::format("%20.10f %20.10f %20.10f\n", sx, sy, sz);
    }
    double pressure = (scs(0, 0) + scs(1, 1) + scs(2, 2)) / 3.0 * unit_transform;
    GlobalV::ofs_running << table << std::endl;
    if (name == "TOTAL-STRESS")
    {
        GlobalV::ofs_running << " TOTAL-PRESSURE: " << std::fixed << std::setprecision(6) << pressure << unit
                             << std::endl
                             << std::endl;
    }
    if (screen)
    {
        std::cout << table << std::endl;
        if (name == "TOTAL-STRESS")
        {
            std::cout << " TOTAL-PRESSURE: " << std::fixed << std::setprecision(6) << pressure << unit << std::endl
                      << std::endl;
        }
    }
    return;
}

}// namespace ModuleIO