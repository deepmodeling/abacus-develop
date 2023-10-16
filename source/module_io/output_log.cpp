#include "output_log.h"
#include "module_base/global_variable.h"
#include "module_base/constants.h"

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
    ModuleBase::GlobalFunc::NEW_PART(name);

    ofs_running << " " << std::setw(8) << "atom" << std::setw(15) << "x" << std::setw(15) << "y" << std::setw(15) << "z"
                << std::endl;
    ofs_running << std::setiosflags(std::ios::showpos);
    ofs_running << std::setprecision(8);

    double fac = 1.0;

    if (!ry)
    {
        fac = ModuleBase::Ry_to_eV / 0.529177;
    }

    if (GlobalV::TEST_FORCE)
    {
        std::cout << " --------------- " << name << " ---------------" << std::endl;
        std::cout << " " << std::setw(8) << "atom" << std::setw(15) << "x" << std::setw(15) << "y" << std::setw(15)
                  << "z" << std::endl;
        std::cout << std::setiosflags(std::ios::showpos);
        std::cout << std::setprecision(6);
    }

    int iat = 0;
    for (int it = 0; it < cell.ntype; it++)
    {
        for (int ia = 0; ia < cell.atoms[it].na; ia++)
        {
            std::stringstream ss;
            ss << cell.atoms[it].label << ia + 1;

            ofs_running << " " << std::setw(8) << ss.str();
            if (std::abs(force(iat, 0)) > output_acc)
                ofs_running << std::setw(20) << force(iat, 0) * fac;
            else
                ofs_running << std::setw(20) << "0";
            if (std::abs(force(iat, 1)) > output_acc)
                ofs_running << std::setw(20) << force(iat, 1) * fac;
            else
                ofs_running << std::setw(20) << "0";
            if (std::abs(force(iat, 2)) > output_acc)
                ofs_running << std::setw(20) << force(iat, 2) * fac;
            else
                ofs_running << std::setw(20) << "0";
            ofs_running << std::endl;

            if (GlobalV::TEST_FORCE)
            {
                std::cout << " " << std::setw(8) << ss.str();
                if (std::abs(force(iat, 0)) > output_acc)
                    std::cout << std::setw(20) << force(iat, 0) * fac;
                else
                    std::cout << std::setw(20) << "0";
                if (std::abs(force(iat, 1)) > output_acc)
                    std::cout << std::setw(20) << force(iat, 1) * fac;
                else
                    std::cout << std::setw(20) << "0";
                if (std::abs(force(iat, 2)) > output_acc)
                    std::cout << std::setw(20) << force(iat, 2) * fac;
                else
                    std::cout << std::setw(20) << "0";
                std::cout << std::endl;
            }

            iat++;
        }
    }

    ofs_running << std::resetiosflags(std::ios::showpos);
    std::cout << std::resetiosflags(std::ios::showpos);
    return;
}

void print_stress(const std::string& name, const ModuleBase::matrix& f, const bool screen, bool ry)
{
    const double output_acc = 1.0e-8;
    GlobalV::ofs_running << " --------------------------- " << name << " ----------------------------" << std::endl;

    double fac = 1.0;

    if (!ry)
    {
        fac = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
    }

    std::cout << std::setprecision(5);
    std::cout << std::setiosflags(std::ios::showpos);

    if (screen)
    {
        std::cout << " ------------------- " << name << " --------------------" << std::endl;
    }

    for (int i = 0; i < 3; i++)
    {
        GlobalV::ofs_running << std::setw(15) << " ";
        if (std::abs(f(i, 0)) > output_acc)
            GlobalV::ofs_running << std::setw(15) << f(i, 0) * fac;
        else
            GlobalV::ofs_running << std::setw(15) << "0";
        if (std::abs(f(i, 1)) > output_acc)
            GlobalV::ofs_running << std::setw(15) << f(i, 1) * fac;
        else
            GlobalV::ofs_running << std::setw(15) << "0";
        if (std::abs(f(i, 2)) > output_acc)
            GlobalV::ofs_running << std::setw(15) << f(i, 2) * fac;
        else
            GlobalV::ofs_running << std::setw(15) << "0";
        GlobalV::ofs_running << std::endl;

        if (screen)
        {
            if (std::abs(f(i, 0)) > output_acc)
                std::cout << std::setw(15) << f(i, 0) * fac;
            else
                std::cout << std::setw(15) << "0";
            if (std::abs(f(i, 1)) > output_acc)
                std::cout << std::setw(15) << f(i, 1) * fac;
            else
                std::cout << std::setw(15) << "0";
            if (std::abs(f(i, 2)) > output_acc)
                std::cout << std::setw(15) << f(i, 2) * fac;
            else
                std::cout << std::setw(15) << "0";
            std::cout << std::endl;
        }
    }

    std::cout << std::resetiosflags(std::ios::showpos);

    return;
}

// print total stress
void printstress_total(const ModuleBase::matrix& scs, bool ry)
{
    double unit_transform = 1;
    if (!ry)
    {
        unit_transform = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
    }

    GlobalV::ofs_running << std::setprecision(6) << std::setiosflags(std::ios::showpos)
                         << std::setiosflags(std::ios::fixed) << std::endl;
    ModuleBase::GlobalFunc::NEW_PART("TOTAL-STRESS (KBAR)");
    std::cout << " ><><><><><><><><><><><><><><><><><><><><><><" << std::endl;
    std::cout << " TOTAL-STRESS (KBAR):" << std::endl;
    std::cout << " ><><><><><><><><><><><><><><><><><><><><><><" << std::endl;

    for (int i = 0; i < 3; i++)
    {
        std::cout << " " << std::setw(15) << scs(i, 0) * unit_transform << std::setw(15) << scs(i, 1) * unit_transform
                  << std::setw(15) << scs(i, 2) * unit_transform << std::endl;

        GlobalV::ofs_running << " " << std::setw(23) << scs(i, 0) * unit_transform << std::setw(23)
                             << scs(i, 1) * unit_transform << std::setw(23) << scs(i, 2) * unit_transform << std::endl;
    }
    double pressure = (scs(0, 0) + scs(1, 1) + scs(2, 2)) / 3.0 * unit_transform;
    std::cout << " TOTAL-PRESSURE: " << pressure << " KBAR" << std::endl;
    GlobalV::ofs_running << " TOTAL-PRESSURE: " << pressure << " KBAR" << std::endl;
    GlobalV::ofs_running << std::setiosflags(std::ios::left);
    std::cout << std::resetiosflags(std::ios::showpos);

    return;
}

}// namespace ModuleIO